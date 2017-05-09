Hetero2 : main.cpp matrix.cpp matrix.h parameters.cpp parameters.h sim_RALRAS.cpp sim_RALRAS.h simulation.cpp simulation.h tool_box.h user_options.h user_options.cpp definitions.h jacobi_eigenvalue.cpp jacobi_eigenvalue.h
	g++ -o Hetero2 main.cpp matrix.cpp parameters.cpp sim_RALRAS.cpp simulation.cpp user_options.cpp jacobi_eigenvalue.cpp
clean:
	rm -f Hetero2
