set(GROFF_COMPILE nroff)
set(GROFF_ROOT ${PROJECT_SOURCE_DIR}/manpages)
ADD_CUSTOM_TARGET (groffdocs ALL)

file(GLOB groffsource ${GROFF_ROOT}/*.nroff)

foreach (gs ${groffsource})
    get_filename_component(mpname ${gs} NAME_WE)
    ADD_CUSTOM_COMMAND(
        SOURCE    ${gs}
        COMMAND   ${GROFF_COMPILE}
        ARGS      ${gs} -c -C -mandoc > ${mpname}.1
        TARGET    groffdocs
        OUTPUTS   ${PROJECT_BINARY_DIR}/manpages/${mpname}.1
    )
    ADD_CUSTOM_COMMAND(
        SOURCE    groffdocs
        TARGET    groffdocs
        DEPENDS   ${PROJECT_BINARY_DIR}/manpages/${mpname}.1
    )
endforeach(gs)

string(REPLACE "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}" gs1 "${groffsource}")
string(REPLACE ".nroff" ".1" gsb "${gs1}")

install(FILES ${gsb} DESTINATION ${MANDIR})
