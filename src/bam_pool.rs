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

        THREAD_READER.with(|cell| {
            let mut reader_opt = cell.borrow_mut();

            // Create reader if this thread doesn't have one yet
            if reader_opt.is_none() {
                *reader_opt = Some(self.create_reader());
            }

            // Use the reader
            f(reader_opt.as_mut().unwrap())
        })
    }
}

// BamReaderPool is Sync because:
// - bam_path and fasta_path are immutable Strings
// - creation_lock is Mutex<()> which is Sync
// Thread-local storage handles the actual readers safely
unsafe impl Sync for BamReaderPool {}
