# Called at build time via cmake -P to regenerate version.py from version.py.in.
# Variables RFC_VERSION_MAJOR/MINOR/PATCH/POSTFIX are passed with -D on the command line.
configure_file("${input}" "${output}" @ONLY)
