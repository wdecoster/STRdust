# Contributing to STRdust

Thanks for your interest in contributing to STRdust! This document covers the
development setup, workflow, and quality standards for the project.

## Getting Started for Contributors

### Prerequisites

- Rust toolchain (install via [rustup](https://rustup.rs/))
- Git

### First-Time Setup

1. Clone the repository:

```bash
git clone https://github.com/wdecoster/STRdust.git
cd STRdust
```

2. Install development tools and git hooks:

```bash
make setup        # Installs rustfmt, clippy, cargo-audit, cargo-outdated
make install-hooks # Installs pre-commit and pre-push hooks
```

The git hooks will automatically run quality checks before commits and pushes, catching issues early.

### Development Workflow

**Quick checks before committing:**

```bash
make pre-commit   # Runs fmt, clippy, and tests
```

**Full CI simulation before pushing:**

```bash
make ci           # Runs fmt-check, clippy, and tests (same as CI)
```

**Other useful commands:**

```bash
make fmt          # Format code
make fmt-check    # Check formatting without modifying files
make clippy       # Run linter
make test         # Run tests
make build        # Build release binary
make build-musl   # Build static MUSL binary
make docs         # Generate and open documentation
make help         # Show all available targets
```

## Testing

STRdust includes comprehensive tests, including specific tests for the `--pathogenic` functionality. Run tests with:

```bash
cargo test
# or
make test
```

Some tests depend on network access (testing the STRchive download functionality) and are skipped by default. They are gated on environment variables — set the variable and run `cargo test` as usual (no `--ignored` needed):

```bash
TEST_PATHOGENIC_NETWORK=1 cargo test   # cache/download behavior
TEST_PATHOGENIC_FULL=1 cargo test      # full --pathogenic workflow (slower)
```

## Code Quality Standards

### Formatting

Code must be formatted with `rustfmt` using the project's configuration (`.rustfmt.toml`):

- Max line width: 100 characters
- Edition: 2024
- Field init shorthand enabled

The pre-commit hook automatically runs formatting checks.

### Linting

The project uses `cargo clippy` for linting with warnings treated as errors. Some clippy warnings are configured to be allowed in `Cargo.toml`:

- `too_many_arguments`: Allowed because bioinformatics functions often require many parameters for configuration

Run clippy with:

```bash
cargo clippy --all-targets --all-features -- -D warnings
# or
make clippy
```

### Security

Security audits run automatically:

```bash
make audit        # Run cargo-audit for vulnerability scanning
make outdated     # Check for outdated dependencies
```

## Dependency Management

This project uses [Dependabot](https://github.com/dependabot) to automatically keep dependencies up to date. Dependabot is configured to:

- Check for Cargo dependency updates weekly on Mondays
- Check for GitHub Actions updates weekly
- Automatically create pull requests for dependency updates
- Group minor and patch updates together for easier review
- Auto-merge patch updates after tests pass
- Require manual review for major version updates

The Dependabot configuration can be found in [`.github/dependabot.yml`](.github/dependabot.yml).

## Continuous Integration

The project uses GitHub Actions for CI/CD:

- **Test workflow**: Runs on all pushes and pull requests
  - Checks formatting (`cargo fmt --check`)
  - Runs clippy with `-D warnings`
  - Runs full test suite
  - Separate job for MUSL static binary build and test
  - Uses cargo caching for faster builds

- **Security workflow**: Runs weekly and on every push/PR
  - Security audit with `cargo-audit`
  - Outdated dependency checks
  - License and security policy enforcement with `cargo-deny`

- **Dependabot workflow**: Automatically tests and merges safe dependency updates

- **Publish workflow**: Creates releases for Linux and macOS when tags are pushed

All CI checks can be simulated locally with `make ci` before pushing.
