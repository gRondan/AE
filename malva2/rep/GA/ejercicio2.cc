#include "newGA.hh"
#include "newGA.req.cc"
#include "newGA.pro.cc"

int main(int argc, char **argv)
{

  char* rutaInstancia = argv[1];
  printf("%s\n",rutaInstancia);
	newGA::Problem pbm;
  pbm.loadInstance(pbm, rutaInstancia);
  newGA::Operator_Pool pool(pbm);
	newGA::SetUpParams cfg(pool);
  // ifstream f1("malva/rep/GA/newGA.cfg");
  ifstream f1("newGA.cfg");
  f1 >> cfg;

	newGA::Solver_Seq solver(pbm,cfg);
	solver.run();

	if (solver.pid()==0)
	{
		solver.show_state();
		cout << "Solution: " << solver.global_best_solution()
		     << " Fitness: " << solver.global_best_solution().fitness() << endl;
		cout << "\n\n :( ---------------------- THE END --------------- :) ";

		ofstream fexit(argv[3]);
		if(!fexit) show_message(13);
		fexit << solver.userstatistics();
    // ofstream myfile;
    // myfile.open (rutaSolucion);
    // std::vector<std::string> empleados;
    // bool flags[pbm.cantidadempleados()];
    // for (int k=1;k<= pbm.cantidadempleados();k++){
    //   string aux = std::to_string(k);
    //   std::string empleado = "e" + aux;
    //   empleados.push_back(empleado);
    //   flags[k-1] = false;
    // }
    // for (int i=0;i<solver.best_solution_trial().array_var().size();i++){//recorro tupla solucion(tareas)
    //   int indice = solver.best_solution_trial().var(i); //Indice de empleado
    //   string aux2 = std::to_string(i+1);
    //   empleados[indice-1] += " t" + aux2;
    //   flags[indice-1] = true;
    // }
    // for (unsigned int i = 0; i < empleados.size(); ++i) {
    //     if(flags[i])
    //       myfile << empleados[i] << "\n";
    // }
    // myfile.close();
	}

	return(0);

}
