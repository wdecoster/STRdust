# Makefile for STRdust development

.PHONY: all build test clean fmt clippy build-musl help

# Default target
all: fmt clippy test build

# Build the project
build:
	cargo build --release

# Build a statically-linked Linux binary using MUSL
build-musl:
	@echo "Building static MUSL binary (x86_64-unknown-linux-musl)"
	rustup target add x86_64-unknown-linux-musl >/dev/null 2>&1 || true
	@if command -v cross >/dev/null 2>&1; then \
		echo "Using cross for reproducible musl build"; \
		OPENSSL_STATIC=1 LIBZ_SYS_STATIC=1 BZIP2_STATIC=1 ZSTD_STATIC=1 LZMA_API_STATIC=1 CURL_STATIC=1 \
		CC_x86_64_unknown_linux_musl=musl-gcc CXX_x86_64_unknown_linux_musl=musl-g++ \
		cross build --release --target x86_64-unknown-linux-musl; \
	else \
		echo "Using cargo. Ensure musl-gcc is available (sudo apt-get install musl-tools)"; \
		OPENSSL_STATIC=1 LIBZ_SYS_STATIC=1 BZIP2_STATIC=1 ZSTD_STATIC=1 LZMA_API_STATIC=1 CURL_STATIC=1 \
		CC_x86_64_unknown_linux_musl=musl-gcc CXX_x86_64_unknown_linux_musl=musl-g++ \
		cargo build --release --target x86_64-unknown-linux-musl; \
	fi
	@echo "Binary: target/x86_64-unknown-linux-musl/release/STRdust"

# Run tests
test:
	cargo test

# Clean build artifacts
clean:
	cargo clean

# Format code
fmt:
	cargo fmt

# Run clippy
clippy:
	cargo clippy --all-targets --all-features -- -D warnings

help:
	@echo "Targets:" \
	&& echo "  build        - Build in release mode" \
	&& echo "  build-musl   - Build static MUSL binary" \
	&& echo "  test         - Run tests" \
	&& echo "  fmt          - Format code" \
	&& echo "  clippy       - Run clippy" \
	&& echo "  clean        - Clean artifacts"
