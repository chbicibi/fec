@echo off

:: このファイルはGA最適化ライブラリのコンパイル用です.
:: Windows環境でmakeコマンドが使用できない場合に代替として使ってください.
:: 2018.7.9_作成

setlocal

:: パラメータ
set lib_name=libfec.a
set target_dir=..\lib

:: コンパイル
call :compile interface
call :compile util
call :compile individual
call :compile moead_unit
call :compile ga_unit
call :compile problem
call :compile basic_optimizer
call :compile moead
call :compile moeadc
call :compile moead_init
call :compile moeadi
call :compile soga
call :compile nsga2
call :compile nsga2c
call :compile tnsdm
call :compile nsga2_init
call :compile nsga2i
call :compile soga_init
call :compile sogai
call :compile ga_unit_selection
call :compile ga_unit_crossover
call :compile ga_unit_mutation
call :compile matrix
call :compile kriging

:: リンク
call :exec "ar cr %lib_name% moeadc.o moead.o moead_init.o moeadi.o moead_unit.o tnsdm.o nsga2c.o nsga2.o nsga2_init.o nsga2i.o soga.o soga_init.o sogai.o ga_unit.o ga_unit_selection.o ga_unit_crossover.o ga_unit_mutation.o basic_optimizer.o individual.o kriging.o problem.o matrix.o util.o interface.o"

:: ライブラリをコピー
call :exec "copy /y %lib_name% %target_dir%"
call :exec "copy /y *.mod %target_dir%"
call :exec "copy /y *.smod %target_dir%"

endlocal
exit /b

:compile
setlocal
set cmd_base=gfortran -O2 -Wall -Wno-unused-dummy-argument -Wno-unused-function -c -o
set cmd=%cmd_base% %~1.o %~1.f08
echo %cmd%
%cmd%
endlocal
exit /b

:exec
echo %~1
%~1
exit /b
