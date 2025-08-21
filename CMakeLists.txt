cmake_minimum_required(VERSION 3.16)
project(uSignal VERSION 0.1.0 LANGUAGES CXX)
enable_testing()

# Some options
option(WRAP_PYTHON "Compile the Python bindings" OFF)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
include(CheckCXXCompilerFlag)

if (${WRAP_PYTHON})
   message("Will build Python bindings")
   find_package(pybind11 REQUIRED)
endif()
find_package(Catch2 3 REQUIRED)
find_package(IPP)
#find_package(FindIPP REQUIRED)
#set(MKL_INTERFACE lp64)
#set(MKL_LINK intel64)
find_package(MKL)

if (${IPP_FOUND})
   message("IPP found")
   add_compile_definitions(WITH_IPP)
else()
   message("Using native filter implementations")
endif()
if (${MKL_FOUND})
   message("MKL found")
   add_compile_definitions(WITH_MKL)
else()
   message("Using native FFT/interpolation implementations")
endif()

##########################################################################################
#                                        Sort Out Python Wrappers                        #
##########################################################################################
#if (WRAP_PYTHON)
#   add_library(pyusignal MODULE
#               python/pyusignal.cpp
##               python/filterRepresentations.cpp
##               python/transforms.cpp
##               python/waveform.cpp
#   )
#   target_link_libraries(pyusignal PRIVATE pybind11::pybind11 pybind11::lto uSignal)
#   set_target_properties(pyusignal PROPERTIES
#                         PREFIX ""
#                         CXX_STANDARD 20
#                         CXX_STANDARD_REQUIRED YES
#                         CXX_EXTENSIONS NO)
#endif()

##########################################################################################
#                                      Define Library Dependencies                       #
##########################################################################################

configure_file(${CMAKE_SOURCE_DIR}/src/private/version.hpp.in
               ${CMAKE_SOURCE_DIR}/include/uSignal/version.hpp)
set(LIBALL MKL::MKL ${IPP_LIBRARIES})
set(LIBALL_PYTHON uSignal ${LIBALL} ${PYTHON_LIBRARIES})

#SET(PYTHON_SRC src/modules/wrap.cpp)
set(UTILITIES_SRCS
    src/system.cpp
    src/vector.cpp
    src/demean.cpp
    src/taper.cpp
    src/window.cpp
    src/utilities/version.cpp
    src/utilities/math/polynomial.cpp
    src/filterRepresentations/finiteImpulseResponse.cpp
    src/filterRepresentations/infiniteImpulseResponse.cpp
    src/filterRepresentations/secondOrderSections.cpp
    src/filterImplementations/infiniteImpulseResponse.cpp)
set(SRCS ${UTILITIES_SRCS})

set(BUILD_SHARED_LIBS YES)
add_library(uSignal SHARED ${SRCS})
target_sources(uSignal
               PUBLIC
               FILE_SET HEADERS
               BASE_DIRS
                  ${CMAKE_CURRENT_SOURCE_DIR}
               FILES
                  include/uSignal/version.hpp
                  include/uSignal/system/system.hpp
                  include/uSignal/filterRepresentations/finiteImpulseResponse.hpp
                  include/uSignal/filterRepresentations/infiniteImpulseResponse.hpp
                  include/uSignal/filterRepresentations/secondOrderSections.hpp
                  include/uSignal/filterImplementations/infiniteImpulseResponse.hpp
                  include/uSignal/demean.hpp
                  include/uSignal/detrend.hpp
                  include/uSignal/taper.hpp
               )
target_include_directories(uSignal
                           PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>)
target_link_libraries(uSignal
                      PRIVATE MKL::MKL ${IPP_LIBRARIES})
set_target_properties(uSignal PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES 
                      CXX_EXTENSIONS NO) 

#########################################################################################
#                                        Unit Tests                                     #
#########################################################################################

add_executable(coreTests
               testing/filterRepresentations/finiteImpulseResponse.cpp
               testing/filterRepresentations/infiniteImpulseResponse.cpp
               testing/filterRepresentations/secondOrderSections.cpp
               testing/utilities/math.cpp
               testing/filterImplementations/infiniteImpulseResponse.cpp
               testing/utilities/demean.cpp
               testing/utilities/window.cpp
               testing/utilities/taper.cpp
)
set_target_properties(coreTests PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES 
                      CXX_EXTENSIONS NO)
target_link_libraries(coreTests PRIVATE uSignal Catch2::Catch2 Catch2::Catch2WithMain)
target_include_directories(coreTests
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)
add_test(NAME coreTests
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing
         COMMAND coreTests)

# Need to copy the answers to the test directory
file(COPY ${CMAKE_SOURCE_DIR}/testing/data DESTINATION .)
#file(COPY ${CMAKE_SOURCE_DIR}/python/benchmarks.py DESTINATION .)

#========================================================================================#
#                                      Installation                                      #
#========================================================================================#
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/${PROJECT_NAME}Config.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION "${version}"
    COMPATIBILITY AnyNewerVersion
)
if (${WRAP_PYTHON})
   install(TARGETS uSignal #pyuSignal
           EXPORT ${PROJECT_NAME}-targets
           FILE_SET HEADERS
           ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
           COMPONENT libraries)
else()
   install(TARGETS uSignal
           EXPORT ${PROJECT_NAME}-targets
           FILE_SET HEADERS 
           ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
           COMPONENT libraries)
endif()
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/uSignal
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
export(EXPORT ${PROJECT_NAME}-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")
install(EXPORT ${PROJECT_NAME}-targets
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
##########################################################################################
#                                     CPack Packaging                                    #
##########################################################################################
include(CPackComponent)
set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_VENDOR "UUSS")
set(CPACK_PACKAGE_CONTACT "ben.baker@utah.edu")
set(CPACK_PACKAGE_LICENSE "MIT")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A seismic signals processing library used for real-time and post-processing applications at University of Utah Seismograph Stations.")
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
set(CPACK_VERBATIM_VARIABLES YES)
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
cpack_add_component(libraries
                    DISPLAY_NAME "uSignal common library")
cpack_add_component(headers
                    DISPLAY_NAME "uSignal library header files"
                    DEPENDS libraries)
set(CPACK_GENERATOR ZIP TGZ)
if (WIN32)
   list(APPEND CPACK_GENERATOR WIX)
elseif (APPLE)
   list(APPEND CPACK_GENERATOR productbuild)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
   find_program(RPMBUILD_PROGRAM rpmbuild)
   if (RPMBUILD_PROGRAM)
      list(APPEND CPACK_GENERATOR RPM)
      #set(CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES abc >= 123")
   endif()
   if (${BUILD_FOR_DEBIAN})
      find_program(DPKG_PROGRAM dpkg REQUIRED)
   else()
      find_program(DPKG_PROGRAM dpkg)
   endif()
   if (DPKG_PROGRAM)
      list(APPEND CPACK_GENERATOR DEB)
      set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}")
      set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON) 
   endif()
else()
   message("Unknown operating system")
endif()
set(CPACK_SOURCE_IGNORE_FILES
  /\\.git/
  \\.swp
  \\.orig
  /CMakeLists\\.txt\\.user
  /private/
)
include(CPack) # Put this last!
