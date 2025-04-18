CFLAGS = -O3 -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -pthread -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

a.out:	main.o algorithm.cpp functions.cpp solution.cpp reduce_sum.cpp residual.cpp
		g++ $^ $(CFLAGS) -o $@

main.o: main.cpp all_includes.h matrix_operations.h parallel_utils.h function_types.h common_types.h
		g++ -c $(CFLAGS) main.cpp

clean:
		rm -f *.out *.o
		