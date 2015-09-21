# - Try to find LibXML2 headers and libraries.
#
# Usage of this module as follows:
#
#     find_package(LibXML2)
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#  LibXML2_ROOT_DIR  Set this variable to the root installation of
#                    LibXML2 if the module has problems finding
#                    the proper installation path.
#                    (it checks also LIBXML2_DIR if variable unset)
#
# Variables defined by this module:
#
#  LIBXML2_FOUND              System has LibXML2 libs/headers
#  LibXML2_LIBRARIES          The LibXML2 libraries
#  LibXML2_LIBRARIES_DIR          The LibXML2 libraries location
#  LibXML2_INCLUDE_DIR        The location of LibXML2 headers

# if no LibXML2_ROOT_DIR is provided try to read an env variable LIBXML2_DIR
IF(NOT DEFINED LibXML2_ROOT_DIR)
    set(libxml2envset $ENV{LIBXML2_DIR})
    IF(libxml2envset)
        set(LibXML2_ROOT_DIR $ENV{LIBXML2_DIR})
    ENDIF(libxml2envset)
ENDIF(NOT DEFINED LibXML2_ROOT_DIR)

find_path(LibXML2_ROOT_DIR
    NAMES include/libxml2/libxml/tree.h
)

find_library(LibXML2_LIBRARIES
    NAMES xml2
    HINTS ${LibXML2_ROOT_DIR}/lib
)

find_path(LibXML2_INCLUDE_DIR
    NAMES libxml/tree.h
    HINTS ${LibXML2_ROOT_DIR}/include/libxml2
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibXML2 DEFAULT_MSG
    LibXML2_LIBRARIES
    LibXML2_INCLUDE_DIR
)

GET_FILENAME_COMPONENT(LibXML2_LIBRARIES_DIR ${LibXML2_LIBRARIES} PATH)

mark_as_advanced(
    LibXML2_ROOT_DIR
    LibXML2_LIBRARIES
    LibXML2_LIBRARIES_DIR
    LibXML2_INCLUDE_DIR
)
