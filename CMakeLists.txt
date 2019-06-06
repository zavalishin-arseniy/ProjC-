# CMakeList.txt: проект CMake для C++Test; включите исходный код и определения,
# укажите здесь логику для конкретного проекта.
#
#cmake_minimum_required (VERSION 3.8)

# Добавьте источник в исполняемый файл этого проекта.
#add_executable (C++Test "C++Test.cpp" )

# TODO: Добавьте тесты и целевые объекты, если это необходимо.


cmake_minimum_required(VERSION 2.8)
project(ProjC-)

set(SOURCE_EXE main.cpp)

set(SOURCE_LIB calc.cpp)

add_library(calc STATIC ${SOURCE_LIB})

add_executable(main ${SOURCE_EXE})

target_link_libraries(main calc)