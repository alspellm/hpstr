
# Declare processing module
module( 
    NAME event 
#    EXTRA_SOURCES EventDict.cxx
    EXTERNAL_DEPENDENCIES ROOT LCIO 
)

#root_generate_dictionary(EventDict ${event_INCLUDE_DIR}/EventDef.h MODULE ${PROJECT_NAME} LINKDEF ${event_INCLUDE_DIR}/EventLinkDef.h)
root_generate_dictionary(EventDict ${event_INCLUDE_DIR}/EventDef.h MODULE ${PROJECT_NAME}  LINKDEF ${event_INCLUDE_DIR}/EventLinkDef.h)

# install ROOT pcm file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libevent_rdict.pcm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libevent.rootmap DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
