#CLion currently cannot diagnose headers correctly if they are not present in the target list for a compilation unit
#So if CLion is running Cmake we add them as additional sources. Normally this is not needed because compilers
#will notice the #include and process them automatically

if($ENV{CLION_IDE})
    file(GLOB HEADERS *.h)
    target_sources(${CURRENT_TARGET} PUBLIC ${HEADERS})
endif()