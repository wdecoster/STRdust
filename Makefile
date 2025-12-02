# Makefile for STRdust development

.PHONY: all build test clean fmt clippy build-musl musl help install-hooks pre-push pre-commit ci setup audit outdated fmt-check docs install bench

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

# Alias for build-musl
musl: build-musl

# Run tests
test:
	cargo test

# Clean build artifacts
clean:
	cargo clean

# Format code
fmt:
	cargo fmt

# Check formatting
fmt-check:
	cargo fmt --check

# Run clippy
clippy:
	cargo clippy --all-targets --all-features -- -D warnings

# Security audit
audit:
	cargo audit

# Check for outdated dependencies
outdated:
	cargo outdated --root-deps-only

# Generate documentation
docs:
	cargo doc --no-deps --open

# Install locally
install:
	cargo install --path .

# Run all checks (CI simulation)
ci: fmt-check clippy test
	@echo "All CI checks passed!"

# Development setup
setup:
	rustup component add rustfmt clippy
	cargo install cargo-audit cargo-outdated

# Install git hooks for automated checks
install-hooks:
	@echo "Installing git hooks..."
	@if [ ! -d .git ]; then \
		echo "❌ Not a git repository (missing .git directory)."; \
		echo "   Run 'git init' or clone the repo with git to enable hooks."; \
		exit 1; \
	fi
	@mkdir -p .git/hooks
	@cp -f .githooks/pre-commit .git/hooks/pre-commit
	@cp -f .githooks/pre-push .git/hooks/pre-push
	@chmod +x .git/hooks/pre-commit .git/hooks/pre-push
	@echo "✅ Git hooks installed successfully!"
	@echo "💡 The hooks will now run automatically on commit and push"

# Run all pre-push checks manually
pre-push: fmt clippy
	@echo "🎉 All pre-push checks passed!"

# Check everything is ready for commit
pre-commit: fmt clippy test
	@echo "Ready for commit!"

# Benchmark (if benchmarks exist)
bench:
	cargo bench

help:
	@echo "Targets:" \
	&& echo "  build        - Build in release mode" \
	&& echo "  build-musl   - Build static MUSL binary" \
	&& echo "  test         - Run tests" \
	&& echo "  fmt          - Format code" \
	&& echo "  fmt-check    - Check formatting without fixing" \
	&& echo "  clippy       - Run clippy linter" \
	&& echo "  audit        - Run security audit" \
	&& echo "  outdated     - Check for outdated dependencies" \
	&& echo "  docs         - Generate and open documentation" \
	&& echo "  clean        - Clean artifacts" \
	&& echo "  ci           - Run all CI checks locally" \
	&& echo "  setup        - Install dev tools" \
	&& echo "  install-hooks - Install git hooks" \
	&& echo "  pre-commit   - Run pre-commit checks" \
	&& echo "  pre-push     - Run pre-push checks"
