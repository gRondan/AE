include ../../environment

all: MainLan MainSeq ejercicio2

clean:
	rm -f MainLan MainSeq MainWan ejercicio2 *.o *% *~

MainLan: newGA.req.o newGA.pro.o MainLan.o
	$(CXX) $(LDFLAGS) $^ $(LOADLIBES) $(CPPFLAGS) -o $@

MainWan: newGA.req.o newGA.pro.o MainWan.o
	$(CXX) $(LDFLAGS) $^ $(LOADLIBES) $(CPPFLAGS) -o $@

MainSeq: newGA.req.o newGA.pro.o MainSeq.o
	$(CXX) $(LDFLAGS) $^ $(LOADLIBES) $(CPPFLAGS) -o $@

LAN:
	$(RUN) -np 10 ./MainLan Config.cfg

WAN:
	$(RUN) -v  -p4pg pgfileWan MainWan

SEQ:
	./MainSeq newGA.cfg ../../ProblemInstances/ONEMAX-instances/onemax10.txt res/sol.txt
