try:
    import pyarrow
except ModuleNotFoundError:
    pyarrow = None


def ensure_pyarrow_compat():
    """
    Awkward versions used by HiggsDNA still look for ``pyarrow.lib.PyExtensionType``
    when reading parquet metadata. Newer pyarrow releases dropped that alias in favor
    of ``ExtensionType``. Re-introduce the old name when needed so parquet I/O keeps
    working without forcing an environment downgrade.
    """
    if pyarrow is None:
        return None

    extension_type = getattr(pyarrow, "ExtensionType", None)
    pyextension_type = getattr(pyarrow, "PyExtensionType", None)

    if pyextension_type is None and extension_type is not None:
        pyarrow.PyExtensionType = extension_type

    pyarrow_lib = getattr(pyarrow, "lib", None)
    if pyarrow_lib is not None:
        lib_pyextension_type = getattr(pyarrow_lib, "PyExtensionType", None)
        lib_extension_type = getattr(pyarrow_lib, "ExtensionType", None)
        if lib_pyextension_type is None and lib_extension_type is not None:
            pyarrow_lib.PyExtensionType = lib_extension_type

    return pyarrow
