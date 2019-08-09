
#ifndef SRC_ELECTROPHYSIOLOGY_IONICMODELS_IONICMODEL_HPP_
#define SRC_ELECTROPHYSIOLOGY_IONICMODELS_IONICMODEL_HPP_


#include <vector>
#include "Util/Factory.hpp"
#include "Electrophysiology/IonicModels/IonicModelsOptions.hpp"


class GetPot;

namespace BeatIt
{

// Only nodal solve
class IonicModel
{

public:

	/// Create a factory
	typedef Factory<IonicModel, std::string>     IonicModelFactory;

	//! Ionic Model base constructor
	/*!
	 *  \param [in] numVar Total number of variables, including the potential
	 *  \param [in] numGatingVar  Number of gating variables (it may be useful for RushLarsen schemes)
	 */
    IonicModel(int numVar, int numGatingVar, const std::string& name = "empty", CellType cell_type = CellType::MCell);
    //! Virtual destructor
    virtual ~IonicModel() {}
    //! Solve method
	/*!
	 *  \param [in] variables Vector containing the local value of all variables (Variables  includes potential)
	 *  \param [in] appliedCurrent value of the applied current
	 *  \param [in] dt        Timestep
	 */
    virtual void solve(std::vector<double>& variables, double appliedCurrent = 0.0, double dt = 1e-3);

	//! Update all the variables in the ionic model
	/*!
	 *  \param [in] variables Vector containing the local value of all variables  (Variables  includes potential)
	 *  \param [in] dt        Timestep
	 */
    virtual void updateVariables(std::vector<double>& variables, double appliedCurrent, double dt) = 0;
    virtual void updateVariables(std::vector<double>& /*v_n*/, std::vector<double>& /*v_np1*/,  double /*appliedCurrent*/, double /*dt*/) {}

    virtual void updateVariables(std::vector<double>& variables, std::vector<double>& /*rhs*/, double appliedCurrent, double dt, bool /*overwrite*/);
    virtual bool isSecondOrderImplemented()
    {
        return false;
    }

	//! Update all the variables in the ionic model
	/*!
     *  \param [in] V transmember potential
	 *  \param [in] variables Vector containing the local value of all variables (Variables does not include potential)
	 *  \param [in] dt        Timestep
	 */
    virtual  void updateVariables(double /*V*/, std::vector<double>& /*variables*/, double /*dt*/) {}

	//! Evaluate total ionic current for the computation of the potential
	/*!
	 *  \param [in] variables Vector containing the local value of all variables (Variables  includes potential)
	 *  \param [in] appliedCurrent value of the applied current
	 *  \param [in] dt        Timestep
	 */
    virtual double evaluateIonicCurrent(std::vector<double>& variables, double appliedCurrent = 0.0, double dt = 0.0) = 0;
    virtual double evaluateIonicCurrentH(std::vector<double>& variables, double appliedCurrent = 0.0, double dt = 0.0, double h = 0.0);

    virtual double evaluatedIonicCurrent(std::vector<double>& /*variables*/, double /*appliedCurrent = 0.0*/, double /*dt = 0.0*/, double /*h = 0.0*/)
    {
        return 0.0;
    }
    virtual double evaluateIonicCurrentTimeDerivative(std::vector<double>& /*variables*/, std::vector<double>& /*old_variables*/, double /*dt = 0.0*/, double /*h = 0.0*/)
    {
        return 0.0;
    }

    virtual double evaluateIonicCurrent(std::vector<double>& /*v_n*/, std::vector<double>& /*rhs*/,
                                        double /*appliedCurrent = 0.0*/, double /*dt = 0.0*/)
    {
        throw std::runtime_error("Calling Base Class IonicModel::evaluateIonicCurrent");
        return 0.0;
    }

    virtual double evaluateSAC(double /*v*/ , double /*I4f*/)
    {
        return 0.0;
    }

	//! Evaluate total ionic current for the computation of the potential
	/*!
     *  \param [in] V transmember potential
	 *  \param [in] variables Vector containing the local value of all variables (Variables  does not  include potential)
	 *  \param [in] appliedCurrent value of the applied current
	 *  \param [in] dt        Timestep
	 */
    virtual double evaluateIonicCurrent(double V, std::vector<double>& variables, double appliedCurrent = 0.0, double dt = 0.0)  = 0;
	//! Initialize the values of the variables
	/*!
	 *  \param [in] variables Vector containing the local value of all variables
	 */
	virtual void initialize(std::vector<double>& variables) = 0;

	//! Initialize the output files with the names
	/*!
	 *  \param [in] output output file where we will save the values
	 */
    virtual void initializeSaveData(std::ostream& output) = 0;

    virtual void setup(GetPot& /*data*/, std::string /*section*/);

    // Number of variables inlcluding transmembrane potential
    int numVariables()
    {
    	return M_numVariables;
    }

    std::string variableName(int numVar)
    {
        if(numVar < M_numVariables-1) return M_variablesNames[numVar];
        else throw std::runtime_error("IonicModel: Accessing name of variable not defined");
    }

    const std::vector<std::string>& variablesNames() const
    {
        return M_variablesNames;
    }

    const std::string& ionicModelName() const
    {
    	return M_ionicModelName;
    }

    virtual void setCellType(CellType type)
    {
        M_cellType = type;
    }

    virtual double membraneCapacitance()
    {
        return M_membrane_capacitance;
    }
    virtual void get_currents(std::vector<double>& currents)
    {

    }
protected:

//    /// Transmembrane potential
//    double M_transmembranePotential;
//    /// Applied current
//    double M_appliedCurrent;
//    /// variables
//    std::vector<double> M_variables;
    /// Total number of variables, including the potential
    int    M_numVariables;
    /// Number of gating variables
    int    M_numGatingVariables;
    /// Name of the variables (excluding the potential)
    std::vector<std::string> M_variablesNames;
    /// typoe of cell
    CellType M_cellType;
    /// Name of the ionic model
    std::string M_ionicModelName;

    /// physical constants:gas constant
    constexpr static double R = 8314.472;
    /// physical constants: Temperature
    constexpr static double T = 310.0;
    /// physical constants: Faraday's constatnt
    constexpr static double F = 96485.3415;
    /// physical constants: R*T/F
    constexpr static double RTONF=(R*T)/F;

    // Overwrite gating variables
    // This is useful for Forward-Backwart Euler
    // and other time integrators
    bool M_overwrite;

    double M_membrane_capacitance;

};



}// BeatIt


#endif /* SRC_ELECTROPHYSIOLOGY_IONICMODELS_IONICMODEL_HPP_ */

