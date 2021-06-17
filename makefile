CXX := g++-10
CXXFLAGS += -Wall -std=c++20 -O3 -fopenmp

LIBSRCS = $(filter-out main.cpp,$(shell find -name \*.cpp))

LIBOBJS = $(patsubst %.cpp,%.o,$(LIBSRCS))

main: main.o libproject.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -L. -o $@ $< -l project

$(LIBOBJS): %.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

libproject.a: $(LIBOBJS)
	ar rcs $@ $^

clean:
	rm -f libproject.a $(LIBOBJS) main

.PHONY: clean
