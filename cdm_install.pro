    TEMPLATE = subdirs
    SUBDIRS += qwt qwt3d cdmlib cdm tests/test1 tests/test2 tests/test3 tests/test4
    CONFIG += ordered

    qwt.file = qwt-6.1.2/qwt.pro
    qwt3d.file = qwtplot3d/qwtplot3d.pro
    cdmlib.file = cdmlib/cdmlib.pro
    cdm.file = cdm/cdm.pro
    tests.file = tests/test1/test1.pro
    tests.file = tests/test2/test2.pro
    tests.file = tests/test3/test3.pro
    tests.file = tests/test4/test4.pro



