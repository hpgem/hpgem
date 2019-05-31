FIND_PACKAGE(Doxygen)
if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
            "Doxygen is needed to build the documentation. Please install it correctly or turn the option")
endif()

if (DOXYGEN_FOUND)

    #This is the configure file for normal doxygen builds
    configure_file(conf/doxygen/doxygen.conf
            ${PROJECT_BINARY_DIR}/conf/doxygen/doxygen.conf  @ONLY IMMEDIATE)

    #The next four and for website doxygen builds. The should be hinded from the public cmake at some point
    configure_file(conf/doxygen/web_doxygen.conf
            ${PROJECT_BINARY_DIR}/conf/doxygen/web_doxygen.conf @ONLY IMMEDIATE)

    configure_file(conf/doxygen/hpg.css
            ${PROJECT_BINARY_DIR}/conf/doxygen/hpg.css @ONLY IMMEDIATE)

    configure_file(conf/doxygen/new_footer.html
            ${PROJECT_BINARY_DIR}/conf/doxygen/new_footer.html @ONLY IMMEDIATE)

    configure_file(conf/doxygen/new_header.html
            ${PROJECT_BINARY_DIR}/conf/doxygen/new_header.html @ONLY IMMEDIATE)





    #-- Add custom targets to both make and clean (delete) the documentation
    add_custom_target 	(doc
            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/conf/doxygen/doxygen.conf
            SOURCES ${PROJECT_BINARY_DIR}/conf/doxygen/doxygen.conf)

    add_custom_target	(docClean
            COMMAND rm -r ${PROJECT_BINARY_DIR}/docs/*
            COMMENT "Cleaning (deleting) the documentation"	)

    add_custom_target	(docWeb
            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/conf/doxygen/web_doxygen.conf
            SOURCES ${PROJECT_BINARY_DIR}/conf/doxygen/web_doxygen.conf)

    add_custom_target	(docPublishAlpha
            COMMAND rsync -vr docs/html/* WebAdmin@einder.ewi.utwente.nl:/home/www/html/hpGEM/assets/doxygen/alpha/
            COMMENT "Publushing (uploading) the curent documentation as ALPHA realese" )
endif()