for %%i in (1.19.3  1.19.5  1.21.6  1.22.4  1.23.5) do (
	cmake -S../.. -B../../build -DRFC_EXPORT_MEX=0 -DRFC_NUMPY_VERSION=%%i -DCMAKE_BUILD_TYPE=Release
	cmake --build ../../build/src/python -t rfcnt --config Release
)
cmake -S../.. -B../../build -DRFC_EXPORT_MEX=0 -DRFC_NUMPY_VERSION= -DCMAKE_BUILD_TYPE=Release
pause
