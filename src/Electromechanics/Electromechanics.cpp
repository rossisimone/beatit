/*
 ============================================================================

 .______    _______     ___   .___________.    __  .___________.
 |   _  \  |   ____|   /   \  |           |   |  | |           |
 |  |_)  | |  |__     /  ^  \ `---|  |----`   |  | `---|  |----`
 |   _  <  |   __|   /  /_\  \    |  |        |  |     |  |     
 |  |_)  | |  |____ /  _____  \   |  |        |  |     |  |     
 |______/  |_______/__/     \__\  |__|        |__|     |__|     
 
 BeatIt - code for cardiovascular simulations
 Copyright (C) 2016 Simone Rossi

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ============================================================================
 */

/**
 * \file Electromechanics.cpp
 *
 * \class Electromechanics
 *
 * \brief This class provides a simple factory implementation
 *
 * For details on how to use it check the test_factory in the testsuite folder
 *
 *
 * \author srossi
 *
 * \version 0.0
 *
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Dec 23, 2016
 *
 */

#include "Electromechanics/Electromechanics.hpp"
//#include "Elasticity/MixedElasticity.hpp"
//#include "Electrophysiology/Monodomain/Monowave.hpp"
#include "Util/IO/io.hpp"
#include "libmesh/numeric_vector.h"
#include "libmesh/exodusII_io.h"
#include "Electrophysiology/IonicModels/IonicModel.hpp"

namespace BeatIt
{

typedef libMesh::TransientExplicitSystem           ActivationSystem;

Electromechanics::Electromechanics( libMesh::EquationSystems& es,
                                    std::string system_name  )
    : M_datafile()
    , M_monowave()
    , M_elasticity()
    , M_equationSystems(es)
    , M_myName(system_name)
    , M_mesh(es.get_mesh())
{
}



Electromechanics::~Electromechanics()
{
}



void
Electromechanics::setup( GetPot& data,
                         std::string electro_section,
                         std::string elasticity_section,
                         std::string activation_section)
{
    std::cout << "* ELECTROMECHANICS: calling setup" << std::endl;

    std::string electro_data = data("electrophysiology", "monodomain.beat");
    std::string elastic_data = data("elasticity", "elasticity.beat");
    std::string  active_data = data("activation", "activation.beat");
    std::cout << "* ELECTROMECHANICS: electrophysiology datafile: " <<  electro_data << std::endl;
    std::cout << "* ELECTROMECHANICS: elastic datafile: " <<  elastic_data << std::endl;
    std::cout << "* ELECTROMECHANICS: activation datafile: " <<  active_data << std::endl;

    // Electrophysiology
    std::cout << "* ELECTROMECHANICS: creating electrophysiology solver " << std::endl;
    GetPot Edata(electro_data);
    M_monowave.reset( new Monowave(M_equationSystems) );
    M_monowave->setup(Edata, electro_section);
    std::string ionic_model_name = M_monowave->M_ionicModelPtr->ionicModelName();

    // Mechanics
    std::cout << "* ELECTROMECHANICS: creating mechanical solver " << std::endl;
    GetPot Mdata(elastic_data);
    M_elasticity.reset( new MixedElasticity(M_equationSystems, "elasticity") );
    M_elasticity->setup(Mdata, elasticity_section);
    // Activation
    std::cout << "* ELECTROMECHANICS: creating activation solver " << std::endl;
    GetPot Adata(active_data);
    M_datafile = Adata;

    //Read output folder from datafile
    std::string output_folder = M_datafile(activation_section+"/output_folder",  "Output");
    M_outputFolder = "./" + output_folder + "/";

    ActivationSystem& activation_system = M_equationSystems.add_system<ActivationSystem>("activation");
    activation_system.add_variable( "gammaf", libMesh::FIRST);
    activation_system.init();

    if(ionic_model_name != "NashPanfilov")
    {
		ActivationSystem& NL06_system = M_equationSystems.add_system<ActivationSystem>("NL06");
		NL06_system.add_variable( "TCa", libMesh::FIRST);
		NL06_system.add_variable( "TCas", libMesh::FIRST);
		NL06_system.add_variable( "Ts", libMesh::FIRST);
		NL06_system.add_variable( "X", libMesh::FIRST);
		NL06_system.add_variable( "F", libMesh::FIRST);
		NL06_system.init();
    }

    BeatIt::createOutputFolder(activation_system.get_mesh().comm(), M_outputFolder);
    M_equationSystems.print_info(std::cout);


    M_exporter.reset( new EXOExporter( M_equationSystems.get_mesh() ) );
    M_exporter->write_equation_systems (M_outputFolder+"em.exo", M_equationSystems);
    M_exporter->append(true);

}

void
Electromechanics::init(double time)
{
    M_monowave->init(time);
}


void
Electromechanics::compute_activation(double dt)
{
//			std::cout << "\n Solving activation: "  <<  std::endl;

    std::string ionic_model_name = M_monowave->M_ionicModelPtr->ionicModelName();
    ActivationSystem& activation_system = M_equationSystems.get_system<ActivationSystem>("activation");
    typedef libMesh::ExplicitSystem                     ParameterSystem;
    ParameterSystem& dummy_system       = M_equationSystems.get_system<ParameterSystem>("dumb");
    typedef libMesh::TransientExplicitSystem           IonicModelSystem;
	IonicModelSystem& ionic_model_system =  M_equationSystems.add_system<IonicModelSystem>("ionic_model");
    typedef libMesh::ExplicitSystem                     ParameterSystem;
    ParameterSystem& I4f_system        = M_equationSystems.get_system<ParameterSystem>("I4f");

    double mu_a = 5.0;

    if(ionic_model_name == "NashPanfilov")
    {
        typedef libMesh::TransientExplicitSystem           WaveSystem;
        WaveSystem& wave_system =  M_equationSystems.get_system<WaveSystem>("wave");

        *dummy_system.solution *= 0.0;
        *dummy_system.solution += *wave_system.solution;
        auto first_local_index = dummy_system.solution->first_local_index();
        auto last_local_index = dummy_system.solution->last_local_index();
        auto i = first_local_index;
         for( ; i < last_local_index; )
         {
             auto val = (*dummy_system.solution)(i);
             if( val < 0.0 )
             {
                 dummy_system.solution->set(i, 0.0);
             }
             i++;

         }
         dummy_system.solution->close();
        *dummy_system.solution *= -0.1;
        *dummy_system.solution -= *activation_system.solution;
        *dummy_system.solution *= dt / mu_a;

        *activation_system.solution += *dummy_system.solution;
    }
    else if (ionic_model_name == "Grandi11")
    {
		ActivationSystem& NL06_system = M_equationSystems.add_system<ActivationSystem>("NL06");

		int num_nl06_vars = 5;


        double Y1 = 39 ; // uM/s
    	double Z1 = 30 ; // /s
		double Y2 = 1.3 ; // /s
		double Z2 = 1.3 ; // /s
		double Y3 = 30 ; // /s
        double Z3 = 1560 ; // uM/s
        double Y4 = 40 ; // /s
        double Yd = 9 ; // s/uM^2
        double Tt = 70 ; // uM
        double B = 1200 ; // /s
        double hc = 0.005 ; // um
        double La = 1.17 ; // um
        double R = 20 ; // /um^2
        double L0 = 0.97 ; // um

        double A = 944.5815 ; // mN/mm^2/um/uM (adjusted)
        A = 0.0625*A ; // mN/mm^2/um/uM (adjusted)
        double B_s = 0.4853 ; // mN/mm^2 (adjusted)
        double h = 3.0843/B_s ; // cm (Kass89)
			auto first_local_index = activation_system.solution->first_local_index();
        auto last_local_index = activation_system.solution->last_local_index();
        auto i = first_local_index;

        int var_num = 35;
        double Cai = 0.0;
        double Cai_rescale = 1e-4;
        double Cai_diastolic = 8.597401e-5/Cai_rescale;
        double I4f = 1.0;
		int num_vars = ionic_model_system.n_vars();
        int var_index = 0;

        double TCa_eff = 0.0;
        double T, TCa, TCas, Ts, X;
        double dTCa, dTCas, dTs, dX;
        double Qb, Qa, Qr, Qd, Qd1, Qd2;
		double F, gf;
         for( ; i < last_local_index; )
         {
//        	         for(int nv = 0; nv < num_vars; ++nv)
//					{
//						var_index =  i * num_vars+ nv;
//						if(nv == 35)
//						{
//						 double cai = (*ionic_model_system.solution)(var_index);
//						 std::cout <<  "nv: " << nv << " = "<< cai << std::endl;
//						}
//					}
//
        	 	 var_index = i*num_vars + var_num;
        	 	 Cai =  (*ionic_model_system.solution)(var_index);
        	 	 Cai /= Cai_rescale;
        	 	 if(Cai <= Cai_diastolic) Cai = 0.0;
        	 	 else Cai -= Cai_diastolic;
				 I4f = (*I4f_system.solution)(i);
				 TCa   = (*NL06_system.solution)(i*num_nl06_vars);
				 TCas = (*NL06_system.solution)(i*num_nl06_vars+1);
				 Ts      = (*NL06_system.solution)(i*num_nl06_vars+2);
				 X       = (*NL06_system.solution)(i*num_nl06_vars+3);
				 F       = (*NL06_system.solution)(i*num_nl06_vars+4);
				 gf = (*activation_system.solution)(i);

				 I4f = 1.0;
				 TCa_eff  = TCa * std::exp(-R*(I4f - La)*(I4f - La)) ;
				 T = Tt - TCa - TCas - Ts;

				 dX = B * (I4f - X - hc); // dX/dt

				 Qb = Y1 * Cai * T - Z1 * TCa;
				 Qa = Y2 * TCa_eff - Z2 * TCas;
				 Qr = Y3 * TCas - Z3 * Ts * Cai ;
				 Qd = Y4 * Ts ;
				 Qd1 = Yd * dX * dX * Ts ;
				 Qd2 = Yd * dX * dX * TCas ;

				 dTCa = Qb - Qa ; // dTCa/dt
				 dTCas = Qa - Qr - Qd2 ; // dTCa*/dt
				 dTs = Qr - Qd - Qd1 ; // dT*/dt


				 TCa +=  dt * dTCa;
				 TCas+= dt * dTCas;
				 Ts += dt * dTs;
				 X += dt * dX;

                F = A * ( TCas + Ts) * (I4f-X);

                gf += dt * ( - F - 10.0 * gf);
//                std::cout << "Cai: " << Cai << ", gf: "  <<  gf << ", TCa: " << TCa << ", TCas: " << TCas << ", Ts: " << Ts << ", X: " << X << ", F: " << F << std::endl;
                activation_system.solution->set(i, gf);
                NL06_system.solution->set(i*num_nl06_vars, TCa);
                NL06_system.solution->set(i*num_nl06_vars+1, TCas);
                NL06_system.solution->set(i*num_nl06_vars+2, Ts);
                NL06_system.solution->set(i*num_nl06_vars+3, X);
                NL06_system.solution->set(i*num_nl06_vars+4, F);
                ++i;
         }
         NL06_system.solution->close();
         activation_system.solution->close();
    }
//    (*activation_system.solution) *= 0.0;
//    (*activation_system.solution) += (*wave_system.solution);
//    (*activation_system.solution) *= -0.14;
//    auto max = activation_system.solution->max();
//    activation_system.solution->add(-max);
}


void
Electromechanics::save_exo(int step, double time)
{
    std::cout << "* ELECTROMECHANICS: EXODUSII::Exporting em.exo at time "   << time << " in: "  << M_outputFolder << " ... " << std::flush;

    M_exporter->write_timestep(  M_outputFolder+"em.exo"
                               , M_equationSystems
                               , step, time );
    std::cout << "done " << std::endl;

}

void
Electromechanics::solve_mechanics()
{
    ActivationSystem& activation_system = M_equationSystems.add_system<ActivationSystem>("activation");
    double dt = 0.0;
    M_elasticity->newton(dt, activation_system.solution.get());
    M_elasticity->evaluate_nodal_I4f();
}

void
Electromechanics::solve_reaction_step( double dt,
                                       double time,
                                       int    step,
                                       bool   useMidpoint,
                                       const  std::string& mass)
{
    typedef libMesh::ExplicitSystem                     ParameterSystem;
    ParameterSystem& I4f_system        = M_equationSystems.get_system<ParameterSystem>("I4f");

    M_monowave->solve_reaction_step(dt, time, step, useMidpoint, mass, I4f_system.solution.get() );
}



} /* namespace BeatIt */
