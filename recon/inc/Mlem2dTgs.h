/**
	\file Mlem2dTgs.h
	\brief Archivo que contiene la definición de la clase Mlem2dTgs. 
	Clase derivada de Mlem, que define el algoritmo Mlem para los sinogramas2d del TGS.

	\todo
	\bug
	\warning
	\author Martín Belzunce (martin.a.belzunce@gmail.com)
	\date 2010.11.11
	\version 1.1.0
*/
#ifndef _MLEM2DTGS_H_
#define _MLEM2DTGS_H_

#include <Mlem.h>
#include <Sinogram2Dtgs.h>
#include <Logger.h>

using namespace::std;


// DLL export/import declaration: visibility of objects
#ifndef LINK_STATIC
	#ifdef WIN32               // Win32 build
		#ifdef DLL_BUILD    // this applies to DLL building
			#define DLLEXPORT __declspec(dllexport)
		#else                   // this applies to DLL clients/users
			#define DLLEXPORT __declspec(dllimport)
		#endif
		#define DLLLOCAL        // not explicitly export-marked objects are local by default on Win32
	#else
		#ifdef HAVE_GCCVISIBILITYPATCH   // GCC 4.x and patched GCC 3.4 under Linux
			#define DLLEXPORT __attribute__ ((visibility("default")))
			#define DLLLOCAL __attribute__ ((visibility("hidden")))
		#else
			#define DLLEXPORT
			#define DLLLOCAL
		#endif
	#endif
#else                         // static linking
	#define DLLEXPORT
	#define DLLLOCAL
#endif



/**
    \brief Clase abstracta del método de reconstrucción MLEM.
    Esta clase abstracta define de forma general una reconstrucción del tipo MLEM. Las clases derivadas
	omplementarán los distintos tipos de reconstrucción, sea 2D, 3D, o con cualquier otra particularidad.Los
	parámetros de reconstrucción son seteados a través de las distintas propiedades de la clase. 
	
    \todo 
*/
#ifdef __cplusplus
//extern "C"
#endif
class DLLEXPORT Mlem2dTgs : public Mlem
{	
	  
	  /// Proyección a reconstruir.
	  /* Objeto del tipo Projection que será la entrada al algoritmo de reconstrucción,
	  puede ser alguno de los distintos tipos de proyección: sinograma 2D, sinograma 3D, etc. */
	  Sinogram2Dtgs* inputProjection;
	  
	  /// Método que calcula la imagen de sensibilidad.
	  /* Método que hace la backprojection de una imagen cosntante para obtener
	  la imagen de sensibilidad necesaria para la reconstrucción. */
	  bool computeSensitivity(Image*);
	  
	public:
		/// Constructores de la clase.
		/* Constructor que carga los parámetros base de una reconstrucción MLEM 2d para el tgs. */
		Mlem2dTgs(Sinogram2Dtgs* cInputProjection, Image* cInitialEstimate, string cPathSalida, string cOutputPrefix, int cNumIterations, int cSaveIterationInterval, bool cSaveIntermediate, bool cSensitivityImageFromFile, Projector* cForwardprojector, Projector* cBackprojector);
		
		/// Constructores de la clase a partir de un archivo de configuración.
		/* Constructor que carga los parámetros base de una reconstrucción MLEM
		a partir de un archivo de configuración con cierto formato dado. */
		Mlem2dTgs(string configFilename);
		
		/// Método que carga los coeficientes de corrección de atenuación desde un archivo interfile para aplicar como corrección.
		/**  Este método habilita la corrección de atenuación y carga la imagen de mapa de atenuación de una imagen interfile.
		      
		*/
		bool setAcfProjection(string acfFilename){return false;};
		
		/// Método que carga un sinograma desde un archivo interfile con la estimación de scatter para aplicar como corrección.
		/**  Este método habilita la corrección por randoms y carga un sinograma para ello.
		*/
		bool setScatterCorrectionProjection(string acfFilename){return false;};
		
		/// Método que carga un sinograma desde un archivo interfile con la estimación de randomc para aplicar como corrección.
		/**  Este método habilita la corrección por randoms y carga un sinograma para ello.
		*/
		bool setRandomCorrectionProjection(string acfFilename){return false;};
		
		/// Método que aplica las correcciones habilitadas según se hayan cargado los sinogramas de atenuación, randoms y/o scatter.
		bool correctInputSinogram(){return false;};
		
		/// Método que realiza la reconstrucción de las proyecciones. 
		bool Reconstruct();
		
};

#endif
