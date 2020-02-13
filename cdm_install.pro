    TEMPLATE = subdirs
    SUBDIRS += qwt qwt3d cdmlib tests cdm
    CONFIG += ordered

    qwt.file = qwt-6.1.2/qwt.pro
    qwt3d.file = qwtplot3d/qwtplot3d.pro
    cdmlib.file = cdmlib/cdmlib.pro
    tests.file = tests/tests.pro
    cdm.file = cdm/cdm.pro

