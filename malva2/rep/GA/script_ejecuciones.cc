#include "malva/rep/GA/newGA.hh"
#include "malva/rep/GA/newGA.req.cc"
#include "malva/rep/GA/newGA.pro.cc"
#include <string.h> 
#include <stdio.h>
#include <stdlib.h>
int main ()
{
	using skeleton newGA;

	system("clear");
        
    int probabilidad_cruzamiento[3];

        probabilidad_cruzamiento[0] = 0.5;
        probabilidad_cruzamiento[1] = 0.75;
        probabilidad_cruzamiento[2] = 1;
	
    int probabilidad_mutacion[3];

        probabilidad_mutacion[0] = 0.001;
        probabilidad_mutacion[1] = 0.01;
        probabilidad_mutacion[2] = 0.1;
    
    int poblacion_total[3];

        poblacion_total[0] = 50;
	    poblacion_total[1] = 125;
        poblacion_total[2] = 200;

    int iteraciones[3];

        iteraciones[0] = 1000;
        iteraciones[1] = 5000;
        iteraciones[2] = 10000;
        
	Problem pbm;
	Operator_Pool pool(pbm);
	SetUpParams cfg(pool);
        
        ofstream archivo_config;
        ofstream archivo_salida;
        archivo_salida.open("malva/rep/GA/salida.txt");
        for(int instancia = 1; instancia <= 5; instancia++){
            string nombre_instancia_string = "malva/rep/GA/instancias/instancia_" + to_string(instancia);
            const char* nombre_instancia_const = nombre_instancia_string.c_str();
            char* nombre_instancia;
            strcpy(nombre_instancia, nombre_instancia_const);
            pbm.loadInstance(nombre_instancia);
            archivo_salida << "--------------INSTANCIA " << to_string(instancia) << "-------------------";
            for(int mutacion = 0; mutacion < 3; mutacion++){
                for(int cruzamiento = 0; cruzamiento < 3; cruzamiento++){
                    for(int poblacion = 0; poblacion < 3; poblacion++){
                        for(int iteracion = 0; iteracion < 3; iteracion++){
                            archivo_salida << "Configuracion: " << "pc: " << cruzamiento << "pm: " << mutacion << "poblacion: " << poblacion;
                            for(int corrida = 1; corrida <= 40; corrida++){
                                archivo_config.open("malva/rep/GA/newGA.cfg");
                                archivo_config << "1    // number of independent runs \n" + to_string(iteracion) + "    // number of generations\n" + to_string(poblacion) + "  // number of individuals\n" + to_string(poblacion) + "  // size of offsprings in each generation\n1 // if replaces parents for offsprings, or only offsprings may be new parents\n1 // display state ?\nSelections	// selections to apply\n1 3   // selection of parents\n2 0 // selection of offsprings\nIntra-Operators // operators to apply in the population\n0 " + to_string(cruzamiento) +"  // crossover & its probability\n1 1.0 " + to_string(mutacion) + "   // mutation & its probability\nInter-Operators // operators to apply between this population and anothers\n0 10 5 1 3 1 5  // operator number, operator rate, number of individuals, selection of indidivual to send and remplace\nLAN-configuration\n1  // refresh global state\n0	// 0: running in asynchronized mode / 1: running in synchronized mode\n1 // interval of generations to check solutions from other populations";
                                archivo_config.close();
                                ifstream f3("malva/rep/GA/newGA.cfg");
                                f3 >> cfg;
                                Solver_Seq solver(pbm,cfg);
                                solver.run();
                                
                                if (solver.pid()==0){
                                    solver.show_state();
                                    archivo_salida << "Corrida: " << corrida;
                                    archivo_salida << " Mejor Fitness: " << solver.global_best_solution().fitness() << endl;
                                    archivo_salida << " Iteracion Mejor Fitness: " << solver.iteration_best_found() << endl;
                                    archivo_salida << " Tiempo para encontrar mejor solucion: " << solver.time_best_found() << endl;
                                    archivo_salida << "**********************************************************";
                                }
                            }
                        }
                    }
                }   
            }
        }
        archivo_salida.close();
      	return(0);

}

