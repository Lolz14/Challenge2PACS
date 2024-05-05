#include <vector>
#include <iostream>
#include "chrono.hpp"
#include "csv.h"
#include "matrix.hpp"
#include <concepts>
#include <complex>


int main()
{

    
    //Type Definition
    using T = double;
    using V = std::complex<double>;

    std::string filename = "lnsp_131.mtx";

    constexpr StorageOrder ro = StorageOrder::RowWise;
    constexpr StorageOrder co = StorageOrder::ColumnWise;

    algebra::Matrix<T,ro> M_RowWise;
    algebra::Matrix<T,co> M_ColWise;

    M_RowWise.reader_mmf(filename);
    M_ColWise.reader_mmf(filename);

    
    std::size_t vector_dimension = M_RowWise.n();

    std::vector<T> vector_;
    vector_.reserve(vector_dimension);

    algebra::Matrix<T,ro> Matrix_V(vector_dimension,1);

    for (size_t i = 0; i < vector_dimension; ++i)
    {
        vector_.emplace_back(i);
        Matrix_V(i,0,i);
    }

    csv output("Output.csv");
    //ROW WISE STORAGE - UNCOMPRESSED
    Timings::Chrono timer1;
    timer1.start();
    M_RowWise*vector_;
    timer1.stop();

    Timings::Chrono timer2;
    timer2.start();
    M_RowWise*Matrix_V;
    timer2.stop();

    Timings::Chrono timer3;
    timer3.start();
    M_RowWise*M_RowWise;
    timer3.stop();

    Timings::Chrono timer4;
    timer4.start();
    M_RowWise.norm<NormType::One>();
    timer4.stop();
    
    Timings::Chrono timer5;
    timer5.start();
    M_RowWise.norm<NormType::Infinity>();
    timer5.stop();

    Timings::Chrono timer6;
    timer6.start();
    M_RowWise.norm<NormType::Frobenius>();
    timer6.stop();
    
    //ROW WISE - COMPRESSED
    M_RowWise.compress();
    Timings::Chrono timer7;
    timer7.start();
    M_RowWise*vector_;
    timer7.stop();
    
    Timings::Chrono timer8;
    timer8.start();
    M_RowWise*Matrix_V;
    timer8.stop();

    Timings::Chrono timer9;
    timer9.start();
    M_RowWise*M_RowWise;
    timer9.stop();

    Timings::Chrono timer10;
    timer10.start();
    M_RowWise.norm<NormType::One>();
    timer10.stop();

    Timings::Chrono timer11;
    timer11.start();
    M_RowWise.norm<NormType::Infinity>();
    timer11.stop();

    Timings::Chrono timer12;
    timer12.start();
    M_RowWise.norm<NormType::Frobenius>();
    timer12.stop();



    //COLWISE - UNCOMPRESSED


    Timings::Chrono timer13;
    timer13.start();
    M_ColWise*vector_;
    timer13.stop();

    Timings::Chrono timer14;
    timer14.start();
    M_ColWise*Matrix_V;
    timer14.stop();

    Timings::Chrono timer15;
    timer15.start();
    M_ColWise*M_ColWise;
    timer15.stop();

    Timings::Chrono timer16;
    timer16.start();
    M_ColWise.norm<NormType::One>();
    timer16.stop();

    Timings::Chrono timer17;
    timer17.start();
    M_ColWise.norm<NormType::Infinity>();
    timer17.stop();

    Timings::Chrono timer18;
    timer18.start();
    M_ColWise.norm<NormType::Frobenius>();
    timer18.stop();


    //COLWISE - COMPRESSED
    M_ColWise.compress();

    Timings::Chrono timer19;
    timer19.start();
    M_ColWise*vector_;
    timer19.stop();

    Timings::Chrono timer20;
    timer20.start();
    M_ColWise*Matrix_V;
    timer20.stop();

    Timings::Chrono timer21;
    timer21.start();
    M_ColWise*M_ColWise;
    timer21.stop();

    Timings::Chrono timer22;
    timer22.start();
    M_ColWise.norm<NormType::One>();
    timer22.stop();

    Timings::Chrono timer23;
    timer23.start();
    M_ColWise.norm<NormType::Infinity>();
    timer23.stop();

    Timings::Chrono timer24;
    timer24.start();
    M_ColWise.norm<NormType::Frobenius>();
    timer24.stop();




    output << "Compressed" <<"Storage"<< "Operation" << "Time (in microseconds)" << row_escape;
 
    output << "False"<<"RowWise"<<"Matrix*vector with std::vector" << timer1.wallTime() << row_escape;
    output << "False"<<"RowWise"<<"Matrix*vector with Matrix<T,O>" << timer2.wallTime() << row_escape;
    output << "False"<<"RowWise"<<"Matrix*Matrix" << timer3.wallTime() << row_escape;
    output << "False"<<"RowWise"<<"Norm-One" << timer4.wallTime() << row_escape;
    output << "False"<<"RowWise"<<"Infinity Norm" << timer5.wallTime() << row_escape;
    output << "False"<<"RowWise"<<"Frobenius Norm" << timer6.wallTime() << row_escape;
    output << "True"<<"RowWise"<<"Matrix*vector with std::vector" << timer7.wallTime() << row_escape;
    output << "True"<<"RowWise"<<"Matrix*vector with Matrix<T,O>" << timer8.wallTime() << row_escape;
    output << "True"<<"RowWise"<<"Matrix*Matrix" << timer9.wallTime() << row_escape;
    output << "True"<<"RowWise"<<"Norm-One" << timer10.wallTime() << row_escape;
    output << "True"<<"RowWise"<<"Infinity Norm" << timer11.wallTime() << row_escape;
    output << "True"<<"RowWise"<<"Frobenius Norm" << timer12.wallTime() << row_escape;
    output << "False"<<"ColWise"<<"Matrix*vector with std::vector" << timer13.wallTime() << row_escape;
    output << "False"<<"ColWise"<<"Matrix*vector with Matrix<T,O>" << timer14.wallTime() << row_escape;
    output << "False"<<"ColWise"<<"Matrix*Matrix" << timer15.wallTime() << row_escape;
    output << "False"<<"ColWise"<<"Norm-One" << timer16.wallTime() << row_escape;
    output << "False"<<"ColWise"<<"Infinity Norm" << timer17.wallTime() << row_escape;
    output << "False"<<"ColWise"<<"Frobenius Norm" << timer18.wallTime() << row_escape;
    output << "True"<<"ColWise"<<"Matrix*vector with std::vector" << timer19.wallTime() << row_escape;
    output << "True"<<"ColWise"<<"Matrix*vector with Matrix<T,O>" << timer20.wallTime() << row_escape;
    output << "True"<<"ColWise"<<"Matrix*Matrix" << timer21.wallTime() << row_escape;
    output << "True"<<"ColWise"<<"Norm-One" << timer22.wallTime() << row_escape;
    output << "True"<<"ColWise"<<"Infinity Norm" << timer23.wallTime() << row_escape;
    output << "True"<<"ColWise"<<"Frobenius Norm" << timer24.wallTime() << row_escape;








  

   return 0;
}

    

   



 


    

 
    

  


    

  
   