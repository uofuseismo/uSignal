from conan import ConanFile
from conan.tools.build import can_run
from conan.tools.cmake import cmake_layout, CMake, CMakeDeps, CMakeToolchain

class uSignalConan(ConanFile):
   name = "uSignal"
   #version = "0.0.1"
   license = "MIT"
   description = "A seismic-focused signals processing library used by various applications at UUSS."
   url = "https://github.com/uofuseismo/uSignal"
   #topics = ("uSignal")
   settings = "os", "compiler", "build_type", "arch"
   options = {"build_tests" : [True, False],}
   default_options = {"hwloc/*:shared": "True",
                      "build_tests" : "True",
                     }
   export_sources = "CMakeLists.txt", "LICENSE", "README.md", "cmake/*", "src/*", "testing/*"
   generators = "CMakeDeps", "CMakeToolchain"

   def requirements(self):
       # dependencies
       self.requires("boost/1.89.0")
       self.requires("onetbb/2022.3.0")

   def build_requirements(self):
       # test dependncies and build tools
       self.test_requires("catch2/3.13.0")

   def layout(self):
       # defines the project layout
       cmake_layout(self)

   def build(self):
       # invokes the build system
       cmake = CMake(self)
       cmake.configure()
       cmake.build()
       #if can_run(self):
       #   # run tests particularly CTest 
       #   cmake.test()

   def test(self):
       if can_run(self):
          cmake.test()

   #def generate(self):
   #    tc = CMakeToolchain(self)
   #    tc.generate()

   def package(self):
       # copies files from the build to package folder
       cmake = CMake(self)
       cmake.install()

   def package_info(self):
       self.cpp_info.libs = ["uSignal"]

