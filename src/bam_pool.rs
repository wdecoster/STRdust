//! Thread-local BAM reader pool for efficient parallel processing
//!
//! This module provides thread-safe BAM reader management for parallel processing.
//! Each thread gets its own BAM reader, avoiding the need for locks during fetch/read
//! operations while ensuring safe reader creation.

use rust_htslib::bam::IndexedReader;
use std::cell::RefCell;
use std::sync::Mutex;

use crate::parse_bam::create_bam_reader;

/// A pool of BAM readers, one per thread
///
/// Uses thread-local storage so each rayon thread gets its own reader.
/// Readers are lazily initialized on first use and reused for subsequent batches.
/// Reader creation is serialized via mutex for htslib safety.
pub struct BamReaderPool {
    bam_path: String,
    fasta_path: String,
    /// Mutex to serialize BAM reader creation (htslib safety requirement)
    creation_lock: Mutex<()>,
}

impl BamReaderPool {
    /// Create a new pool with the given BAM and FASTA paths
    pub fn new(bam_path: String, fasta_path: String) -> Self {
        Self { bam_path, fasta_path, creation_lock: Mutex::new(()) }
    }

    /// Create a BAM reader (serialized for htslib safety)
    fn create_reader(&self) -> IndexedReader {
        // Serialize reader creation to prevent concurrent htslib index operations
        let _guard = self.creation_lock.lock().unwrap();
        create_bam_reader(&self.bam_path, &self.fasta_path)
    }

    /// Execute a closure with a BAM reader for the current thread
    ///
    /// The reader is obtained from thread-local storage (or created if first use).
    /// This is the main entry point for using the pool in parallel code.
    pub fn with_reader<F, T>(&self, f: F) -> T
    where
        F: FnOnce(&mut IndexedReader) -> T,
    {
        // Use thread-local storage for the reader
        // Each rayon thread will have its own reader
        thread_local! {
            static THREAD_READER: RefCell<Option<IndexedReader>> = const { RefCell::new(None) };
        }

        // Take the reader out of the cell and drop the borrow before running the closure.
        //
        // Holding `borrow_mut()` across `f` panics the moment anything re-enters this on the
        // same thread, and rayon makes that reachable: the caller here is already inside a
        // `par_iter`, and `phase_insertions` starts a nested one to build its distance
        // matrix. Work-stealing can schedule that inner task onto the very thread sitting in
        // this borrow, so `--unphased --threads N` panicked with "RefCell already borrowed".
        // It never showed up on phased input because the clustering module is not reached.
        //
        // Taking ownership costs nothing in the common case - the reader is moved back on
        // the way out and reused exactly as before. A re-entrant call now finds the cell
        // empty and builds its own reader rather than panicking, which is one extra reader
        // on a nested path, not one per locus.
        let mut reader = THREAD_READER
            .with(|cell| cell.borrow_mut().take())
            .unwrap_or_else(|| self.create_reader());

        let result = f(&mut reader);

        THREAD_READER.with(|cell| {
            *cell.borrow_mut() = Some(reader);
        });
        result
    }
}

// BamReaderPool is Sync because:
// - bam_path and fasta_path are immutable Strings
// - creation_lock is Mutex<()> which is Sync
// - the readers themselves never cross threads: each lives in thread-local storage, is
//   taken out and put back by the same thread, and is never handed to another
// Note the last point is what `with_reader` has to preserve; an implementation that held a
// reference across a nested call would break re-entrancy, not soundness.
unsafe impl Sync for BamReaderPool {}
