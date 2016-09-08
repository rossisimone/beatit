/*
 * main_factory.cpp
 *
 *  Created on: 13/feb/2015
 *      Author: srossi
 */

#include <memory>
#include <functional>

#include "Util/Factory.hpp"
#include "Util/IO/io.hpp"

/*!
 *    This test shows how to use the Factory class
 *
 *    Assume you have an abstract base class Base defining the interfaces
 *    of the class. From our base class we derive several classes, such as
 *    Derived1, Derived2, Derived3 and so on.
 *    Since all the interfaces are defined in the base class we need to define
 *    only ( a pointer to ) the Base class object,
 *
 *    On the other hand, since the base class is abstract we cannot
 *    create. Therefore, when we instantiate the pointer we will need to create
 *    one of the derived classes.
 *    For example:
 *    std::shared_ptr<Base> BasePtr ( new Derived1 );
 *
 *    Let's assume we want to decide the type at run-time, say Derived7.
 *    We could implement a method crateBasePtr(int ID) such that we can
 *    have a switch function such as the following
 *
 *    switch( ID )
 *    {
 *         case 1:
 *            return new Derived1;
 *         case 2:
 *            return new Derived2;
 *         ....
 *
 *         case 7;
 *            return new Derived7;
 *         default:
 *          exit(-1);
 *    }
 *
 *
 *    Each time we implement a new Derived class we need to update this method.
 *
 *    Using the Factory class we avoid this.
 *
 *    In this example I will show how to actually use the factory class in AthenaVMS
 *    The key point is defining a typedef to the factory
 *    in the base class. Then, for any derived class we need to define:
 *    1) a function creating the object;
 *    2) register the object into the factory (details below).
 *
 *    Note that the class Base, Derived1, Derived2 should be in different files!
 *    I will point out the start of each file in the following.
 *
 *
 *
 *
 */

/*!
 *    FILE: Base.hpp
 *
 *    We start by defining the Base class defining the interfaces
 */

class Base
{
public:

    /*!
     *   This is the key Element to use the Factory
     *
     *   The template parameters refer to:
     *   1) the Base class defining the interfaces.
     *   2) a key for which we can select the different Derived classes
     *   In this example we will use strings to select
     *   which class to use, and this class it will be
     *   of type Base (actually it will be a child of it)
     *
     *   UPDATE:
     *   The factory now takes three template parameters
     *   where the last can be void
     *   This is useful for passing a simple object to the constructor
     *
     */
    typedef Factory<Base, std::string> factory_Type;

    // Abstract interface
    virtual void doSomething() = 0;

    virtual ~Base() {}
};



/*!
 *    FILE: Derived1.hpp
 *
 *    Simple implementation on the class Derived1
 */


class Derived1 : public Base
{
public:
    void doSomething()
    {
        std::cout << "Class Derived1 is lazy! I do not do anything ... " << std::endl;
    };
};

/*!
 *    We need to create a function that creates an object of type Derived1
 *    but it returns a pointer to the Base class
 *    This function will be registered in the factory,
 *    such that when we will create the object,
 *    we will actually be calling this function
 */
Base* createDerived1()
{
    return new Derived1;
}

/*!
 *   The registration is done in an anonymous namespace
 *   With this definition we can create a new object
 *   of type Derived1 using the string "Derived1"
 *   We will need to call the function:
 *   Base::factory_Type::Create( "Derived" );
 */
namespace
{
    static bool registerDerived1 = Base::factory_Type::Register( "Derived1", &createDerived1 );
}

/*!
 *    FILE: Derived2.hpp
 *
 *    Simple implementation on the class Derived2
 *    For details on the implementation,
 *    see the comments on class Derived1 (above)
 */


class Derived2 : public Base
{
public:

    void doSomething()
    {
        std::cout << "Class Derived2 is happy! It's doing something ... " << std::endl;
    }
};


Base* createDerived2()
{
    return new Derived2;
}


namespace
{
    static bool registerDerived2 = Base::factory_Type::Register( "Derived2", &createDerived2 );
}


/*!
 *    FILE: main.cpp ( or where you actually need to use the class )
 *
 */

//IN THIS TEST WE PASS THE Factory
int main(int argc, char *argv[])
{
	BeatIt::printBanner(std::cout);

    std::cout << "Factory test: this test should fail as we are creating an object that does not exists!" << std::endl;

    /// Here we actually create the objects
    std::unique_ptr<Base> obj1;
    obj1.reset( Base::factory_Type::Create("Derived1") );
    std::shared_ptr<Base> obj2 ( Base::factory_Type::Create("Derived2") );

    /// We use the same interface but we have different behavior
    obj1->doSomething();
    obj2->doSomething();

    /// Here we actually create an object that we did not registered. It should bomb!


    try
    {
        Base * obj3 = Base::factory_Type::Create("Derived3");
    }
    catch (const std::runtime_error& error)
    {
    	std::cout << "Caught runtime error: " << std::flush;
    	std::cout << error.what() << std::endl;
        return EXIT_SUCCESS;

    }

    return EXIT_FAILURE;
}
