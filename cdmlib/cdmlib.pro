######################################################################
# Automatically generated by qmake (2.01a) Mon Dec 27 12:07:07 2010
######################################################################

TEMPLATE 	= 	lib

VERSION         =       0.3.5

TARGET 		=       cdmlib

CONFIG          +=      staticlib warn_on 

DEPENDPATH 	+= .

DESTDIR      	= lib

QT      	+=

DEFINES 	+=      CDMVERSION=\\\"$$VERSION\\\"  DEBUG 
CONFIG(fftw) {
DEFINES 	+=      USE_FFTW
}

DEFINES 	+= 	QT_NO_DEBUG_OUTPUT

QMAKE_CC        =       gfortran 

QMAKE_CFLAGS    += -Warray-bounds -fcray-pointer -w -cpp 

QMAKE_LFLAGS    = 

QMAKE_CFLAGS_RELEASE    = -O3 

QMAKE_CFLAGS_THREAD =

HEADERS 	+=      cdmlib.h

SOURCES		+= 	cdmlib.f \
                        calculdate.f \
                        bicgstarplus.f \
                        cors.f \
                        beampropagation.f \
                        beampropagationmacro.f \
                        calculgreenfft.f \
                        calculgreenquadfft.f \
                        calculforcefft.f \
                        calculforcefftopt.f \
			coeffmie.f \
                        comparaisonxyz.f \
                        derivechamp.f \
			d07hre.f \
			d09hre.f \
			d113re.f \
			d132re.f \
			d1mach.f \
			d9lgmc.f \
			dadhre.f \
			dasyjy.f \
			dbesj.f  \
                        dcsevl.f \
			dchhre.f \
			dcuhre.f \
                        derivative.f \
			derivativefield2.f \
			dfshre.f \
			dgamlm.f \
			dgamma.f \
			dinhre.f \
			djairy.f \
			dlngam.f \
			drlhre.f \
			dtrhre.f \
			dznrm2.f \
			espace_libre_derive_mat.f \
			espace_libre.f \
                        fdump.f \
                        fouriertoimage.f \
                        fouriertoimage2.f \                        
                        fftsingletonz.f \                        
			gaussiand.f \
                        gaussian.f \
                        gaussianfft.f \
                        gaussianpara.f \
			gausskronrodpattersonmulti.f \
			gpbicg2.f \
			gpbicgar2.f \
			gpbicgar.f \
			gpbicg.f \
                        gpbicgsafe.f \
                        gpbicgplus.f \                        
                        gpbicor.f \
			i1mach.f \
                        incidentarbitrary.f \
			initds.f \
			interpdielec5.f \
			inversemat33.f \
                        inversemat33r.f \
                        inverserig.f \
                        inverserigopt.f \
			irradiance.f \
			j4save.f \
			local-macro.f \
			module_mie.f \
                        objectarbitrary.f \
			objectcube.f \
			objectcylindre.f \
			objectellipse.f \
			objectnspheres.f \
                        objectpara.f \
                        objectparanxnynz.f \
			objectsphereconcentric.f \
                        objectsphere.f \
                        objectsphererandom.f \
                        objectsphererandomnxnynz.f \                        
                        objectinhomo.f \
                        objectparainhomo.f \
                        objectparainhomonxnynz.f \
                        ondecirced.f \
                        ondecirce.f \
                        ondecircekxky.f \
			ondeplaned.f \
                        ondeplane.f \
                        ondeplanekxky.f \
                        ondeplanemulti.f \
			ondeplanedmulti.f \
			pimzbicgstab.f \
			polaepstens.f \
			polarisabilite.f \
                        produitfftmatvect3.f \
                        produitfftmatvectopt.f \
			qmrbicgstab.f \
			qmrpim.f \
			tfqmr.f \
			wofz.f \
			xercnt.f \
			xerhlt.f \
			xermsg.f \
			xerprn.f \
			xersve.f \
			xgetua.f \
			zcg.f \
			zgedid.f \
                        relecture.f \
                        comparaison.f \
                        diffractefft2d.f \
                        diffractefft2denergie.f \
                        diffractefft2dinc.f \
                        diffractefft2dlens.f \
			computegcfft2d.f \
                        propaesplibreint.f \
                        dipoleinc.f \
                        dipoleincder.f \
                        random.f \
                        microsbf.f \
                        microsdf.f \
                        microssansinc.f \
                        deltakroutine.f \
                        polamodifie.f \
                        hdf5close.f \ 
                        hdf5open.f \
                        hdf5read1d_int.f \
                        hdf5read2d_int.f \
                        hdf5write1d_int.f \
                        hdf5write2d_int.f \
                        hdf5create.f \
                        hdf5read1d.f \
                        hdf5read2d.f \
                        hdf5write1d.f \
                        hdf5write2d.f \
                        writehdf5mic.f

                        
INCLUDEPATH 	+= .

CONFIG(fftw) {
	LIBS 		+= 	-lgfortran -lfftw3_omp -lfftw3 -lm 
} else {
	LIBS 		+= 	-lgfortran -lm 
}

CONFIG(hdf5) {
# sur centos, fedora, etc...
  exists( /usr/lib64/gfortran/modules ) {
	LIBS 		+= 	-I/usr/lib64/gfortran/modules -I/usr/include -L/usr/lib64 -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
  }
# sur ubuntu
  exists( /usr/include/hdf5/serial ) {
	LIBS 		+= 	-I/usr/include/hdf5/serial -I/usr/include -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
  }
} else {
	LIBS 		+= 	
}

QMAKE_DISTCLEAN += lib/*



