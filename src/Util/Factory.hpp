/**
 * \class Factory
 *
 *
 * \brief This class provides a simple factory implementation
 *
 * For details on how to use it check the test_factory in the testsuite folder
 *
 *
 * \author Simone Rossi
 *
 * \version 0.0
 *
 *
 * Contact: bv@bnl.gov
 *
 * Created on: Aug 6, 2016

 *
 */

#ifndef SRC_UTIL_FACTORY_HPP_
#define SRC_UTIL_FACTORY_HPP_

#include <map>
#include <iostream>


template< class Product, class Identifier, class Argument>
class FactoryArg
{
public:
    typedef Product* (*Creator)(Argument &);
private:
    typedef std::map<Identifier, Creator> FactoryMap;
public:

    /*!
     * 	\brief Registers the creator function in the factory and returns 'true' if registration was successful
     *
     *  \param [in] id string/value corresponding to the constructor
     *  \param [in] creator pointer to a function calling the object constructor
     */
    static bool Register(const Identifier& id, Creator creator)
    {
        return getMap().insert( typename FactoryMap::value_type(id, creator) ).second;
    }


    /*!
     * 	\brief Creates the object spcified by the Identifier if
     *
     *  \param [in] id string/value corresponding to the constructor
     */
    static Product* Create(const Identifier& id, Argument& arg)
    {
    	// Check if the id was registered in the factory
        typename FactoryMap::const_iterator i =  getMap().find(id);
    	// If it was registered then call the constructor
        if (i != getMap().end())
        {
            return (i->second)(arg);
        }
        // Else show the available options and throw an error
        else
        {
        	// show available keys
            std::cout << "Factory: Product of type '" << id << "' has not been registered in the factory!!!" << std::endl;
            std::cout << "Available keys: " << std::endl;
            for( auto it = getMap().begin(); it != getMap().end(); ++it)
            {
                std::cout << it->first << std::endl;
            }

            // Throw a runtime error
            throw std::runtime_error("Factory: stopping execution!");

        }
    }

private:
    /*!
     * 	\brief return the static map contained in the factory
     */
    static FactoryMap& getMap()
    {
        static FactoryMap S_map;
        return S_map;
    }

};




template< class Product, class Identifier>
class Factory
{
public:
    typedef Product* (*Creator)();
private:
    typedef std::map<Identifier, Creator> FactoryMap;
public:

    /*!
     *  \brief Registers the creator function in the factory and returns 'true' if registration was successful
     *
     *  \param [in] id string/value corresponding to the constructor
     *  \param [in] creator pointer to a function calling the object constructor
     */
    static bool Register(const Identifier& id, Creator creator)
    {
        return getMap().insert( typename FactoryMap::value_type(id, creator) ).second;
    }


    /*!
     *  \brief Creates the object spcified by the Identifier if
     *
     *  \param [in] id string/value corresponding to the constructor
     */
    static Product* Create(const Identifier& id)
    {
        // Check if the id was registered in the factory
        typename FactoryMap::const_iterator i =  getMap().find(id);
        // If it was registered then call the constructor
        if (i != getMap().end())
        {
            return (i->second)();
        }
        // Else show the available options and throw an error
        else
        {
            // show available keys
            std::cout << "Factory: Product of type '" << id << "' has not been registered in the factory!!!" << std::endl;
            std::cout << "Available keys: " << std::endl;
            for( auto it = getMap().begin(); it != getMap().end(); ++it)
            {
                std::cout << it->first << std::endl;
            }

            // Throw a runtime error
            throw std::runtime_error("Factory: stopping execution!");

        }
    }

private:
    /*!
     *  \brief return the static map contained in the factory
     */
    static FactoryMap& getMap()
    {
        static FactoryMap S_map;
        return S_map;
    }

};



#endif /* SRC_UTIL_FACTORY_HPP_ */
