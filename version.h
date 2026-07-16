// The single source of the FSMP-kMC version. The engine bakes it into
// --version, the Studio parses it for the window title and the start page
// (fsmp_gui.__version__), tools/make_bundle.py names local bundles with it
// and the release workflow refuses a tag that does not match it.
// Bump it before tagging a release (tag v0.5.0 <-> "0.5.0" here).
#define FSMP_VERSION "0.5.0-dev"
