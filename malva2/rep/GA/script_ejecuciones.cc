#include "newGA.hh"
#include "newGA.req.cc"
#include "newGA.pro.cc"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
int main ()
{
	using skeleton newGA;

	// system("clear");

    float probabilidad_cruzamiento[3];

        probabilidad_cruzamiento[0] = 0.5;
        probabilidad_cruzamiento[1] = 0.75;
        probabilidad_cruzamiento[2] = 1;

    float probabilidad_mutacion[3];

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
	SetUpParams cfg2(pool);

  ofstream archivo_config;
  // ofstream archivo_salida;
	ofstream csv;
  // archivo_salida.open("salida.txt");
	csv.open("salida.csv",  fstream::app);

	ifstream f3("newGA.cfg");
	f3 >> cfg;
	f3.close();
  for(int instancia = 1; instancia <= 5; instancia++){
			const char* nombre_instancia_const = "instancias/instancia_";
			char integer_string[32];
			sprintf(integer_string, "%d", instancia);
			size_t len = strlen(nombre_instancia_const);

		  char* ret = new char[len+2];

		  strcpy(ret, nombre_instancia_const);
		  ret[len] = integer_string[0];
		  ret[len+1] = '\0';
			pbm.loadInstance(pbm, ret);
      // for(int mutacion = 0; mutacion < 1; mutacion++){
      //     for(int cruzamiento = 0; cruzamiento < 1; cruzamiento++){
      //         for(int poblacion = 0; poblacion < 1; poblacion++){
								 // for(int iteracion = 0; iteracion < 40; iteracion++){
							   //      archivo_salida << "Configuracion: " << "pc: " << cruzamiento << "pm: " << mutacion << "poblacion: " << poblacion;
								 //      archivo_config.open("newGA.cfg");
									// 		archivo_config << "1    // number of independent runs \n" +
									// 		to_string(100) + "    // number of generations\n" +
									// 		to_string(poblacion_total[poblacion]) + "  // number of individuals\n" +
									// 		to_string(100) +	"  // size of offsprings in each generation\n"+
									// 		"1 // if replaces parents for offsprings, or only offsprings may be new parents\n"+
									// 		"0 // display state ?\n"+
									// 		"Selections	// selections to apply\n"+
									// 		"1 3   // selection of parents\n"+
									// 		"2 0 // selection of offsprings\n"+
									// 		"Intra-Operators // operators to apply in the population\n"+
									// 		"0 " + to_string(probabilidad_cruzamiento[cruzamiento]) +"  // crossover & its probability\n"+
									// 		"1 1.0 " + to_string(probabilidad_mutacion[mutacion]) + "   // mutation & its probability\n"+
									// 		"Inter-Operators // operators to apply between this population and anothers\n"+
									// 		"0 10 5 1 3 1 5  // operator number, operator rate, number of individuals, selection of indidivual to send and remplace\n"+
									// 		"LAN-configuration\n"+
									// 		"1  // refresh global state\n"+
									// 		"0	// 0: running in asynchronized mode / 1: running in synchronized mode\n"+
									// 		"1 // interval of generations to check solutions from other populations";
								 //      archivo_config.close();

											for(int corrida = 1; corrida <= 40; corrida++){
												Solver_Seq solver(pbm,cfg);
									      solver.run();

									      if (solver.pid()==0){
                            solver.show_state();
                            // archivo_salida << "Corrida: " << corrida;
														csv << instancia << ";" <<solver.global_best_solution().fitness() << ";" << solver.iteration_best_found() << ";" << solver.time_best_found() << endl;
                            // archivo_salida << " Mejor Fitness: " << solver.global_best_solution().fitness() << endl;
                            // archivo_salida << " Iteracion Mejor Fitness: " << solver.iteration_best_found() << endl;
                            // archivo_salida << " Tiempo para encontrar mejor solucion: " << solver.time_best_found() << endl;
                            // archivo_salida << "**********************************************************";
                        }
												// Solver_Seq solver2(pbm,cfg2);
									      // solver2.run();
												//
									      // if (solver2.pid()==0){
                        //     solver2.show_state();
                        //     // archivo_salida << "Corrida: " << corrida;
												// 		csv << instancia << ";" <<solver2.global_best_solution().fitness() << ";" << solver2.iteration_best_found() << ";" << solver2.time_best_found() << endl;
                        //     // archivo_salida << " Mejor Fitness: " << solver.global_best_solution().fitness() << endl;
                        //     // archivo_salida << " Iteracion Mejor Fitness: " << solver.iteration_best_found() << endl;
                        //     // archivo_salida << " Tiempo para encontrar mejor solucion: " << solver.time_best_found() << endl;
                        //     // archivo_salida << "**********************************************************";
                        // }
											}

                // }
      //         }
      //     }
      // }
  }
  // archivo_salida.close();
	csv.close();
	return(0);

}
