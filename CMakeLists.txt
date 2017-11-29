#  CMAKE Main File for the AR-PET Image Reconstruction library (APIRL)
#
#  Martin Belzunce, UTN-FRBA, Proyecto AR-PET (CNEA)
#  Copyright (c) 2009

PROJECT(APIRL)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Type of build (Debug or Release). CMAKE_BUILD_TYPE es una variable del tipo cache, por lo que no conviene
# forzarla a un valor, mejor definirla cuando uno llama a cmake.
#set(CMAKE_BUILD_TYPE	Debug)

# Di estoy en modo debug definio una variable para el preprocesador, para poder incorporar
# pedazos de código de debuggeo con compilación condicional:
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-D__DEBUG__)
endif()

# Disable secure Warnings for stdc functions.
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

###################### OPENMP CONFIGURATION #############################
INCLUDE(FindOpenMP)
IF(OPENMP_FOUND)
    OPTION(USE_OPENMP    "Manages the frame grabbing on windows" ON)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	message(STATUS "Flags: ${OpenMP_C_FLAGS} ${OpenMP_CXX_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
     SET(ENV{OMP_STACKSIZE} 512M)
ELSE(OPENMP_FOUND)
    SET(USE_OPENMP OFF)
ENDIF(OPENMP_FOUND)
#######################################################################

################## CAMBIO EL STACK SIZE POR EL OPENMP #########
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --main-stacksize=8500000")
SET (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} --main-stacksize=8500000")
#########################################################

# Genero una variable con los directorios de invlude para cada subproyecto. Las creo en este, para que
# puedan ser accesibles desde cada uno de ellos.
SET(cmd_Headers_Dir ${CMAKE_SOURCE_DIR}/cmd/inc)
SET(recon_Headers_Dir ${CMAKE_SOURCE_DIR}/recon/inc)
SET(utils_Headers_Dir ${CMAKE_SOURCE_DIR}/utils/inc)
SET(data_Headers_Dir ${CMAKE_SOURCE_DIR}/data/inc)
SET(reconGPU_Headers_Dir ${CMAKE_SOURCE_DIR}/reconGPU/inc)

# El directorio de instalación para modo debug es sobre el mismo directorio de build, en un directorio
# llamado instal. Para la versión release no se cambia, se deja el por default, que suele depender
# del sistema operativo:
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	SET(APIRL_INSTALL_PATH "${CMAKE_BINARY_DIR}/debug")
else()
	SET(APIRL_INSTALL_PATH "${CMAKE_INSTALL_PREFIX}")
endif()	
MESSAGE(STATUS "DIRECTORIO DE INSTALACION: ${APIRL_INSTALL_PATH}")



# Variables de entorno necesarias en linux, si no están definidas las fuerzo a un valor por default:
if(CMAKE_HOST_UNIX)
  if(NOT ENV{CUDA_SDK_INSTALL_PATH})
    SET(ENV{CUDA_SDK_INSTALL_PATH} "$ENV{HOME}/NVIDIA_GPU_Computing_SDK")
  endif()
endif()

option(ENABLE_GPU "build CUDA" OFF)
if(ENABLE_GPU)
	# Seteo en ON CUDA_VERBOSE_BUILD para visualizar todos los resultados del nvcc.
	SET(CUDA_VERBOSE_BUILD ON)
	SET(CUDA_BUILD_CUBIN ON)
	# Defino una variable que me indice si las librería para cuda son dinámicas o estáticas. Esto lo hago
	# así porque he tenido algún problema para hacerlas dinámicas.
	SET(CUDA_LIBRARIES_BUILD	STATIC)
	#SET(CUDA_LIBRARIES_BUILD	SHARED)
	# Algunas opciones adicionales para cuda.
	#SET(CUDA_USE_FAST_MATH ON)
	# Tipo de Cache que se habilita:
	#SET(CUDA_USE_CACHE_L1_L2 ON)
	SET(CUDA_USE_CACHE_ONLY_L2 ON)
	#SET(CUDA_MAX_REG_COUNT 32)
	#set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -rdc=true -lineinfo")
endif(ENABLE_GPU)

#set(CUDA_NVCC_FLAGS_RELEASE ${CUDA_NVCC_FLAGS_RELEASE} ";-03")
set(CUDA_NVCC_FLAGS_DEBUG "${CUDA_NVCC_FLAGS_DEBUG} -G")
# Seteo las compute capabilities con las que quiero compilar (Por default es la más baja: 1.0):
set(CUDA_COMPUTE_CAPABILITY 3.5 CACHE STRING "1.0-3,2.0-1,3.0,3.5")
if(CUDA_COMPUTE_CAPABILITY STREQUAL 1.0)
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_10")
elseif(CUDA_COMPUTE_CAPABILITY STREQUAL 1.1)
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_11")
elseif(CUDA_COMPUTE_CAPABILITY STREQUAL 1.2)
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_12")
elseif(CUDA_COMPUTE_CAPABILITY STREQUAL 1.3)
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_13")
elseif(CUDA_COMPUTE_CAPABILITY STREQUAL 2.0)
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_20")
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode=arch=compute_20,code=sm_20")
elseif(CUDA_COMPUTE_CAPABILITY STREQUAL 2.1)
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_21")
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode=arch=compute_21,code=sm_21")
elseif(CUDA_COMPUTE_CAPABILITY STREQUAL 3.0)
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_30")
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode=arch=compute_30,code=sm_30")
elseif(CUDA_COMPUTE_CAPABILITY STREQUAL 3.5)
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_35")
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode=arch=compute_35,code=sm_35")
endif()
# Agregado Martin Belzunce (14/04/11): Flags según el compute capability por default 1.0 (el mínimo):
if(CUDA_USE_FAST_MATH)
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -use_fast_math")
endif()
  if(CUDA_MAX_REG_COUNT)
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -maxrregcount=${CUDA_MAX_REG_COUNT}")
  endif()
if(CUDA_USE_CACHE_L1_L2)
 set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xptxas=-v,-dlcm=ca -lineinfo")
  #set(nvcc_flags ${nvcc_flags} "-Xptxas --dlcm=ca")
elseif(CUDA_USE_CACHE_ONLY_L2)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xptxas=-v,-dlcm=cg -lineinfo")
  #set(nvcc_flags ${nvcc_flags} "-Xptxas --dlcm=cg")
else()
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xptxas=-v -lineinfo")
endif()
MESSAGE(STATUS "NVCC FLAGS: ${CUDA_NVCC_FLAGS}")

################## CHECK IF CUDA AVAILABLE ############################
# Directorio donde se encuentra el FindCUDA:
if(ENABLE_GPU)
	# In cmake 3.0 there is a findcuda embedded in cmake, for previous version use our findcuda.
	if(CMAKE_MAJOR_VERSION LESS 3)
		set(CUDA_NVCC_FLAGS "") # Our findcuda set the flags based on the variables on the top. To avoid redefinition clean it.
		SET(FindCUDA_Dir ${CMAKE_SOURCE_DIR}/FindCUDA)
		set(CMAKE_MODULE_PATH "${FindCUDA_Dir}/CMake/cuda" ${CMAKE_MODULE_PATH})
	endif()
	find_package(CUDA)
	if (CUDA_FOUND)
		SET(ENABLE_GPU ON)
		add_definitions(-D__USE_CUDA__)
	else (CUDA_FOUND)
		SET(ENABLE_GPU OFF)	
	endif (CUDA_FOUND)
endif(ENABLE_GPU)
#######################################################################

# Agrego los directorios de cada proyecto:
add_subdirectory(recon)
add_subdirectory(data)
add_subdirectory(utils)
add_subdirectory(cmd)
if(ENABLE_GPU)
	add_subdirectory(reconGPU)
	add_subdirectory(cmdGPU)
endif(ENABLE_GPU)

# Configuraciones Adicionales de CUDA. En windows compilo en 32 bits al menos que obtenga un
# Visual Studio de 64 bits. En Linux si directamente compilo según el sistema operativo.
# Cuando la variable CUDA_64_BIT_DEVICE_CODE está en ON se compila con la opción -m64 para
# 64 bits, de lo contrario se utiliza la opción -m32.
if(CMAKE_HOST_WIN32 OR CMAKE_HOST_WIN64)
	if(CMAKE_SIZEOF_VOID_P EQUAL 8)
	  set(CUDA_64_BIT_DEVICE_CODE OFF)
	else(CMAKE_SIZEOF_VOID_P EQUAL 4)
	  set(CUDA_64_BIT_DEVICE_CODE OFF)
	endif(CMAKE_SIZEOF_VOID_P EQUAL 8)
else()
	if(CMAKE_SIZEOF_VOID_P EQUAL 8)
	  set(CUDA_64_BIT_DEVICE_CODE ON)
	else(CMAKE_SIZEOF_VOID_P EQUAL 4)
	  set(CUDA_64_BIT_DEVICE_CODE OFF)
	endif(CMAKE_SIZEOF_VOID_P EQUAL 8)
	
endif()

# Check if Doxygen is available
#find_package(Doxygen REQUIRED)
find_package(Doxygen)

OPTION(INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF)
OPTION(USE_DOT "Set to ON to perform diagram generation with graphviz" OFF)
OPTION(USE_LATEX "Set to ON to build latex documentation" OFF)
OPTION(USE_CHM "Set to ON to build CHM Windows documentation" OFF)
SET(INSTALL_DOC ON)
IF (INSTALL_DOC)
	INCLUDE("${CMAKE_SOURCE_DIR}/Doxygen/generateDOC.cmake" )
	GENERATE_DOCUMENTATION(${CMAKE_SOURCE_DIR}/${PROJECT_NAME}.dox.in)
ENDIF()

# Instalo los archivos de sample.
#FILE(GLOB Sample_Files RELATIVE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/samples/*.txt | ${CMAKE_SOURCE_DIR}/samples/*.par)
FILE(GLOB Sample_Files ABSOLUTE ${CMAKE_SOURCE_DIR}/samples/*.txt | ${CMAKE_SOURCE_DIR}/samples/*.par)

Install(FILES ${Sample_Files} DESTINATION "${APIRL_INSTALL_PATH}/Samples")




