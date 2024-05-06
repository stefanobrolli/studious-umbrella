Main: Main.o Basis.o Props.o SelfEn.o Functions.o Basis.h Props.h SelfEn.h Functions.h 
	g++ -O3 -o RichModSCGF  Main.o Basis.o Props.o SelfEn.o Functions.o
Basis.o: Basis.cpp Basis.h
	g++ -O3 -c -o Basis.o Basis.cpp 
Props.o: Props.cpp Props.h
	g++ -O3 -c -o Props.o Props.cpp 
SelfEn.o: SelfEn.cpp SelfEn.h
	g++ -O3 -c -o SelfEn.o SelfEn.cpp 
Functions.o: Functions.cpp Functions.h 
	g++ -O3 -c -o Functions.o Functions.cpp 
Main.o: Main.cpp
	g++ -O3 -c -o Main.o Main.cpp 
clean: 
	rm *.o RichModSCGF
