
add_executable(FemSC main.cpp ${FemSCSOURCE})

add_executable(TestIntegrate TestIntegrate.cpp ${FemSCSOURCE})

add_executable(TestOneDProblem TestOneDProblem.cpp ${FemSCSOURCE})

add_executable(TestTwoDProblem TestTwoDProblem.cpp ${FemSCSOURCE})

add_executable(TestMatrix Testmatrix.cpp ${FemSCSOURCE})

add_executable(Fem2018 main2018.cpp ${FemSOURCE2018} ${SOURCES2018})

