#include "newGA.hh"
#include "newGA.req.cc"
#include "newGA.pro.cc"
#include <string.h> 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
int main ()
{
	using skeleton newGA;

	system("clear");
         
        
	Problem pbm;
	Operator_Pool pool(pbm);
	SetUpParams cfg(pool); 

        ofstream archivo_config;
        ofstream archivo_salida;
        archivo_salida.open("salida.txt");    

        for(int instancia = 1; instancia <= 3; instancia++){ 
            const char* nombre_instancia_const = "instancias/instancia_";
            char integer_string[32];
            sprintf(integer_string, "%d", instancia);
            printf("%s\n", integer_string);
            printf("%s\n", nombre_instancia_const);
            size_t len = strlen(nombre_instancia_const);
            char* ret = new char[len+2];
            strcpy(ret, nombre_instancia_const);
            ret[len] = integer_string[0];
            ret[len+1] = '\0';
            printf("%s\n", ret);
            pbm.loadInstance(pbm, ret);
            archivo_salida << "--------------INSTANCIA " << to_string(instancia) << "-------------------\n";
            ifstream f3("newGA.cfg");
            f3 >> cfg;
            clock_t tStart = clock();
            Solver_Seq solver(pbm,cfg);
            solver.run();
            if (solver.pid()==0){
               solver.show_state();
               archivo_salida << " Fitness: " << solver.global_best_solution().fitness() << endl;
               archivo_salida << " Tiempo para encontrar la solucion: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
               archivo_salida << "**********************************************************\n";
            }       
        }   

        archivo_salida.close();
      	return(0);
}

