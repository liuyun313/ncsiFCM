ncsiFCM: main.o FCM.o
	g++  main.o FCM.o -o ncsiFCM

main.o: main.cpp 
	g++ -c main.cpp -o main.o

FCM.o: FCM.cpp FCM.h
	g++ -c FCM.cpp -o FCM.o

clean:
	rm -rf *.o ncsiFCM

     
