#ifndef JETCORRECTOR_H
#define JETCORRECTOR_H

// SingleJetCorrector
// v3.0
// Author: Yi Chen
// 
// This class applies JEC for any given level using TF1 as the workhorse
// Supposedly runs faster than v1.0
// v3.0: one can add list of text files to apply them one by one
// v4.0: JetUncertainty (v1.0, Yi Chen) merged in, and JES variations:
//       GetCorrectedPT(Variation::UP / DOWN) = corrected pT * (1 +- JEU),
//       GetUncertainty() = {down, up} fractions at the corrected pT
// v4.1: CorrectedPt(...) / Uncertainty(...) take the jet as arguments instead
//       of setters, like JetSmearer and JetSelector

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include <stdexcept>
#include <utility>

// shared with JetSmearer.h
#ifndef JET_VARIATION_ENUM
#define JET_VARIATION_ENUM
enum class Variation { NOMINAL = 0, DOWN = 1, UP = 2 };
#endif

// JetUncertainty
// v1.0
// Author: Yi Chen
// 
// This class gives you jet uncertainties
//

class JetUncertainty
{
private:
   enum Type { TypeNone, TypeJetPT, TypeJetEta, TypeJetPhi, TypeJetArea, TypeRho };
   bool Initialized;
   double JetPT, JetEta, JetPhi, JetArea, Rho;
   std::vector<std::vector<Type>> BinTypes;
   std::vector<std::vector<double>> BinRanges;
   std::vector<std::vector<double>> PTBins;
   std::vector<std::vector<double>> ErrorLow;
   std::vector<std::vector<double>> ErrorHigh;
public:
   JetUncertainty()                  { Initialized = false; }
   JetUncertainty(std::string File)  { Initialized = false; Initialize(File); }
   ~JetUncertainty()                 {}
   void SetJetPT(double value)     { JetPT = value; }
   void SetJetEta(double value)    { JetEta = value; }
   void SetJetPhi(double value)    { JetPhi = value; }
   void SetJetArea(double value)   { JetArea = value; }
   void SetRho(double value)       { Rho = value; }
   void Initialize(std::string FileName);
   std::vector<std::string> BreakIntoParts(std::string Line);
   bool CheckDefinition(std::string Line);
   std::string StripBracket(std::string Line);
   JetUncertainty::Type ToType(std::string Line);
   std::pair<double, double> GetUncertainty();
   double GetValue(Type T);
};

void JetUncertainty::Initialize(std::string FileName)
{
   int nvar = 0, npar = 0;
   std::string CurrentFormula = "";
   std::vector<Type> CurrentDependencies;
   std::vector<Type> CurrentBinTypes;

   std::ifstream in(FileName.c_str());

   while(in)
   {
      char ch[1048576];
      ch[0] = '\0';
      in.getline(ch, 1048577, '\n');

      if(ch[0] == '\0')
         continue;

      bool IsDefinition = CheckDefinition(ch);
      std::vector<std::string> Parts = BreakIntoParts(StripBracket(ch));

      if(Parts.size() == 0)
         continue;

      if(IsDefinition == true)
      {
         // Found a definition line - update current formula

         nvar = atoi(Parts[0].c_str());
         if(Parts.size() <= nvar + 1)
            continue;
         npar = atoi(Parts[nvar+1].c_str());
         if(Parts.size() <= nvar + 1 + npar + 1)
            continue;

         CurrentFormula = Parts[nvar+1+npar+1];
         
         CurrentBinTypes.clear();
         for(int i = 0; i < nvar; i++)
            CurrentBinTypes.push_back(ToType(Parts[1+i]));

         CurrentDependencies.clear();
         for(int i = 0; i < npar; i++)
            CurrentDependencies.push_back(ToType(Parts[nvar+2+i]));
      }
      else
      {
         // Otherwise it's a line with actual JECs, add it to the list

         if(Parts.size() < nvar * 2 + 1)
            continue;

         BinTypes.push_back(CurrentBinTypes);

         std::vector<double> Ranges;
         for(int i = 0; i < nvar * 2; i++)
            Ranges.push_back(atof(Parts[i].c_str()));
         for(int i = 0; i + 1 < (int)Ranges.size(); i = i + 2)
            if(Ranges[i] > Ranges[i+1])
               std::swap(Ranges[i], Ranges[i+1]);
         BinRanges.push_back(Ranges);

         int N = atoi(Parts[nvar*2].c_str());

         std::vector<double> pt, errorlow, errorhigh;

         for(int i = 0; i + 2 < N; i = i + 3)
         {
            pt.push_back(atof(Parts[nvar*2+i+1].c_str()));
            errorlow.push_back(atof(Parts[nvar*2+i+2].c_str()));
            errorhigh.push_back(atof(Parts[nvar*2+i+3].c_str()));
         }
         PTBins.push_back(pt);
         ErrorLow.push_back(errorlow);
         ErrorHigh.push_back(errorhigh);
      }
   }

   in.close();

   Initialized = true;
}

std::vector<std::string> JetUncertainty::BreakIntoParts(std::string Line)
{
   std::stringstream str(Line);

   std::vector<std::string> Result;

   while(str)
   {
      std::string Temp = "";
      str >> Temp;
      if(Temp == "")
         continue;
      Result.push_back(Temp);
   }

   return Result;
}

bool JetUncertainty::CheckDefinition(std::string Line)
{
   for(int i = 0; i < (int)Line.size(); i++)
   {
      if(Line[i] == '{')
         return true;
      if(Line[i] == ' ')
         continue;
      return false;
   }

   return false;
}

std::string JetUncertainty::StripBracket(std::string Line)
{
   for(int i = 0; i < (int)Line.size(); i++)
   {
      if(Line[i] == '{' || Line[i] == '}')
      {
         Line.erase(Line.begin() + i);
         i = i - 1;
      }
   }

   return Line;
}

JetUncertainty::Type JetUncertainty::ToType(std::string Line)
{
   if(Line == "JetPt")    return TypeJetPT;
   if(Line == "JetEta")   return TypeJetEta;
   if(Line == "JetPhi")   return TypeJetPhi;
   if(Line == "JetA")     return TypeJetArea;
   if(Line == "Rho")      return TypeRho;

   std::cerr << "[JetUncertainty] Warning: variable type " << Line << " not found!" << std::endl;

   return TypeNone;
}

std::pair<double, double> JetUncertainty::GetUncertainty()
{
   if(Initialized == false)
      return std::pair<double, double>(-1, -1);

   int N = BinTypes.size();

   for(int iE = 0; iE < N; iE++)
   {
      bool InBin = true;

      for(int iB = 0; iB < (int)BinTypes[iE].size(); iB++)
      {
         double Value = GetValue(BinTypes[iE][iB]);
         if(Value < BinRanges[iE][iB*2] || Value > BinRanges[iE][iB*2+1])
            InBin = false;
      }

      if(InBin == false)
         continue;

      if(PTBins[iE].size() == 0)
         return std::pair<double, double>(-1, -1);

      double JetPT = GetValue(TypeJetPT);

      if(JetPT < PTBins[iE][0])
         return std::pair<double, double>(ErrorLow[iE][0], ErrorHigh[iE][0]);
      if(JetPT >= PTBins[iE][PTBins[iE].size()-1])
         return std::pair<double, double>(ErrorLow[iE][PTBins[iE].size()-1], ErrorHigh[iE][PTBins[iE].size()-1]);
      
      int Bin = 0;
      for(int i = 0; i < (int)PTBins[iE].size() - 1; i++)
      {
         if(JetPT >= PTBins[iE][i] && JetPT < PTBins[iE][i+1])
         {
            Bin = i;
            break;
         }
      }

      double Low = ErrorLow[iE][Bin] + (ErrorLow[iE][Bin+1] - ErrorLow[iE][Bin]) / (PTBins[iE][Bin+1] - PTBins[iE][Bin]) * (JetPT - PTBins[iE][Bin]);
      double High = ErrorHigh[iE][Bin] + (ErrorHigh[iE][Bin+1] - ErrorHigh[iE][Bin]) / (PTBins[iE][Bin+1] - PTBins[iE][Bin]) * (JetPT - PTBins[iE][Bin]);

      return std::pair<double, double>(Low, High);
   }

   return std::pair<double, double>(-1, -1);
}

double JetUncertainty::GetValue(Type T)
{
   if(T == TypeNone)      return 0;
   if(T == TypeJetPT)     return JetPT;
   if(T == TypeJetEta)    return JetEta;
   if(T == TypeJetPhi)    return JetPhi;
   if(T == TypeJetArea)   return JetArea;
   if(T == TypeRho)       return Rho;

   return -1;
}



class JetCorrector;
class SingleJetCorrector;

class SingleJetCorrector
{
private:
   enum Type { TypeNone, TypeJetPT, TypeJetEta, TypeJetPhi, TypeJetArea, TypeRho };
   bool Initialized;
   bool IsFunction;
   double JetPT, JetEta, JetPhi, JetArea, Rho;
   std::vector<std::string> Formulas;
   std::vector<std::vector<double>> Parameters;
   std::vector<std::vector<Type>> BinTypes;
   std::vector<std::vector<double>> BinRanges;
   std::vector<std::vector<Type>> Dependencies;
   std::vector<std::vector<double>> DependencyRanges;
   std::vector<TF1 *> Functions;
public:
   SingleJetCorrector()                  { Initialized = false; }
   SingleJetCorrector(std::string File)  { Initialized = false; Initialize(File); }
   ~SingleJetCorrector()                 { for(auto P : Functions) if(P != nullptr) delete P; }
   void SetJetPT(double value)     { JetPT = value; }
   void SetJetEta(double value)    { JetEta = value; }
   void SetJetPhi(double value)    { JetPhi = value; }
   void SetJetArea(double value)   { JetArea = value; }
   void SetRho(double value)       { Rho = value; }
   void Initialize(std::string FileName);
   std::vector<std::string> BreakIntoParts(std::string Line);
   bool CheckDefinition(std::string Line);
   std::string StripBracket(std::string Line);
   SingleJetCorrector::Type ToType(std::string Line);
   double GetCorrection();
   double GetCorrectedPT();
   double GetValue(Type T);
private:
   std::string Hack4(std::string Formula, char V, int N);
};

class JetCorrector
{
private:
   std::vector<SingleJetCorrector> JEC;
   double JetPT, JetEta, JetPhi, JetArea, Rho;
   JetUncertainty JEU;
   bool HasUncertainty = false;
public:
   JetCorrector()                               {}
   JetCorrector(std::string File)               { Initialize(File); }
   JetCorrector(std::vector<std::string> Files) { Initialize(Files); }
   JetCorrector(std::vector<std::string> Files, std::string UncertaintyFile)
                                                { Initialize(Files); InitializeUncertainty(UncertaintyFile); }
   void Initialize(std::string File)            { std::vector<std::string> X; X.push_back(File); Initialize(X); }
   void Initialize(std::vector<std::string> Files);
   void InitializeUncertainty(std::string File) { JEU.Initialize(File); HasUncertainty = true; }
   void SetJetPT(double value)     { JetPT = value; }
   void SetJetEta(double value)    { JetEta = value; }
   void SetJetPhi(double value)    { JetPhi = value; }
   void SetJetArea(double value)   { JetArea = value; }
   void SetRho(double value)       { Rho = value; }
   double GetCorrection();
   double GetCorrectedPT();
   double GetCorrectedPT(Variation V);
   std::pair<double, double> GetUncertainty();
   double CorrectedPt(double RawPT, double Eta, double Phi, double Rho, double Area,
                      Variation V = Variation::NOMINAL);
   std::pair<double, double> Uncertainty(double RawPT, double Eta, double Phi, double Rho,
                                         double Area);
};

void JetCorrector::Initialize(std::vector<std::string> Files)
{
   JEC.clear();
   for(auto File : Files)
      JEC.push_back(SingleJetCorrector(File));
}

double JetCorrector::GetCorrection()
{
   double PT = GetCorrectedPT();
   if(PT < 0)
      return -1;
   return PT / JetPT;
}

double JetCorrector::GetCorrectedPT()
{
   double PT = JetPT;

   for(int i = 0; i < (int)JEC.size(); i++)  
   {
      JEC[i].SetJetPT(PT);
      JEC[i].SetJetEta(JetEta);
      JEC[i].SetJetPhi(JetPhi);
      JEC[i].SetRho(Rho);
      JEC[i].SetJetArea(JetArea);

      PT = JEC[i].GetCorrectedPT();

      if(PT < 0)
         break;
   }

   return PT;
}

// JES uncertainty {down, up} as fractions, taken at the corrected pT, e.g. to
// shift a smeared pT; {-1, -1} where the uncertainty file has no entry
std::pair<double, double> JetCorrector::GetUncertainty()
{
   if(HasUncertainty == false)
      throw std::runtime_error("JetCorrector: no uncertainty file for GetUncertainty / Variation::UP/DOWN");

   double PT = GetCorrectedPT();
   if(PT < 0)
      return std::pair<double, double>(-1, -1);

   JEU.SetJetPT(PT);
   JEU.SetJetEta(JetEta);
   JEU.SetJetPhi(JetPhi);
   JEU.SetJetArea(JetArea);
   JEU.SetRho(Rho);
   return JEU.GetUncertainty();
}

// argument style: the jet in one call, no setters (sets them, then as above)
double JetCorrector::CorrectedPt(double RawPT, double Eta, double Phi, double Rho, double Area,
                                 Variation V)
{
   SetJetPT(RawPT);
   SetJetEta(Eta);
   SetJetPhi(Phi);
   SetRho(Rho);
   SetJetArea(Area);
   return GetCorrectedPT(V);
}

std::pair<double, double> JetCorrector::Uncertainty(double RawPT, double Eta, double Phi,
                                                    double Rho, double Area)
{
   SetJetPT(RawPT);
   SetJetEta(Eta);
   SetJetPhi(Phi);
   SetRho(Rho);
   SetJetArea(Area);
   return GetUncertainty();
}

// JES variation: corrected pT * (1 + up) or (1 - down); -1 where the
// uncertainty file has no entry
double JetCorrector::GetCorrectedPT(Variation V)
{
   double PT = GetCorrectedPT();
   if(V == Variation::NOMINAL || PT < 0)
      return PT;

   std::pair<double, double> U = GetUncertainty();
   if(V == Variation::UP)
      return (U.second < 0) ? -1 : PT * (1 + U.second);
   return (U.first < 0) ? -1 : PT * (1 - U.first);
}

void SingleJetCorrector::Initialize(std::string FileName)
{
   int nvar = 0, npar = 0;
   std::string CurrentFormula = "";
   std::vector<Type> CurrentDependencies;
   std::vector<Type> CurrentBinTypes;

   std::ifstream in(FileName.c_str());

   while(in)
   {
      char ch[1048576];
      ch[0] = '\0';
      in.getline(ch, 1048577, '\n');

      if(ch[0] == '\0')
         continue;

      bool IsDefinition = CheckDefinition(ch);
      std::vector<std::string> Parts = BreakIntoParts(StripBracket(ch));

      if(Parts.size() == 0)
         continue;

      if(IsDefinition == true)
      {
         // Found a definition line - update current formula

         nvar = atoi(Parts[0].c_str());
         if(Parts.size() <= nvar + 1)
            continue;
         npar = atoi(Parts[nvar+1].c_str());
         if(Parts.size() <= nvar + 1 + npar + 1)
            continue;

         CurrentFormula = Parts[nvar+1+npar+1];
         
         CurrentBinTypes.clear();
         for(int i = 0; i < nvar; i++)
            CurrentBinTypes.push_back(ToType(Parts[1+i]));

         CurrentDependencies.clear();
         for(int i = 0; i < npar; i++)
            CurrentDependencies.push_back(ToType(Parts[nvar+2+i]));
      }
      else
      {
         // Otherwise it's a line with actual JECs, add it to the list

         if(Parts.size() < nvar * 2 + npar * 2 + 1)
            continue;

         std::vector<double> Parameter;
         for(int i = nvar * 2 + npar * 2 + 1; i < (int)Parts.size(); i++)
            Parameter.push_back(atof(Parts[i].c_str()));
         Parameters.push_back(Parameter);

         Dependencies.push_back(CurrentDependencies);

         if(CurrentDependencies.size() == 4)
            Formulas.push_back(Hack4(CurrentFormula, 't', Parameter.size()));
         else
            Formulas.push_back(CurrentFormula);

         std::vector<double> Ranges;
         for(int i = nvar * 2 + 1; i < nvar * 2 + 1 + npar * 2; i++)
            Ranges.push_back(atof(Parts[i].c_str()));
         for(int i = 0; i + 1 < (int)Ranges.size(); i = i + 2)
            if(Ranges[i] > Ranges[i+1])
               std::swap(Ranges[i], Ranges[i+1]);
         DependencyRanges.push_back(Ranges);

         BinTypes.push_back(CurrentBinTypes);

         Ranges.clear();
         for(int i = 0; i < nvar * 2; i++)
            Ranges.push_back(atof(Parts[i].c_str()));
         for(int i = 0; i + 1 < (int)Ranges.size(); i = i + 2)
            if(Ranges[i] > Ranges[i+1])
               std::swap(Ranges[i], Ranges[i+1]);
         BinRanges.push_back(Ranges);

         Functions.push_back(nullptr);
      }
   }

   in.close();

   Initialized = true;
}

std::vector<std::string> SingleJetCorrector::BreakIntoParts(std::string Line)
{
   std::stringstream str(Line);

   std::vector<std::string> Result;

   while(str)
   {
      std::string Temp = "";
      str >> Temp;
      if(Temp == "")
         continue;
      Result.push_back(Temp);
   }

   return Result;
}

bool SingleJetCorrector::CheckDefinition(std::string Line)
{
   for(int i = 0; i < (int)Line.size(); i++)
   {
      if(Line[i] == '{')
         return true;
      if(Line[i] == ' ')
         continue;
      return false;
   }

   return false;
}

std::string SingleJetCorrector::StripBracket(std::string Line)
{
   for(int i = 0; i < (int)Line.size(); i++)
   {
      if(Line[i] == '{' || Line[i] == '}')
      {
         Line.erase(Line.begin() + i);
         i = i - 1;
      }
   }

   return Line;
}

SingleJetCorrector::Type SingleJetCorrector::ToType(std::string Line)
{
   if(Line == "JetPt")    return TypeJetPT;
   if(Line == "JetEta")   return TypeJetEta;
   if(Line == "JetPhi")   return TypeJetPhi;
   if(Line == "JetA")     return TypeJetArea;
   if(Line == "Rho")      return TypeRho;

   std::cerr << "[SingleJetCorrector] Warning: variable type " << Line << " not found!" << std::endl;

   return TypeNone;
}

double SingleJetCorrector::GetCorrection()
{
   if(Initialized == false)
      return -1;

   int N = Formulas.size();

   for(int iE = 0; iE < N; iE++)
   {
      bool InBin = true;

      for(int iB = 0; iB < (int)BinTypes[iE].size(); iB++)
      {
         double Value = GetValue(BinTypes[iE][iB]);
         if(Value < BinRanges[iE][iB*2] || Value > BinRanges[iE][iB*2+1])
            InBin = false;
      }

      if(InBin == false)
         continue;

      if(Dependencies[iE].size() == 0)
         return -1;   // huh?
      if(Dependencies[iE].size() > 4)
      {
         std::cerr << "[SingleJetCorrector] There are " << Dependencies[iE].size() << " parameters!" << std::endl;
         return -1;   // huh?
      }

      double V[3] = {0, 0, 0};
      for(int i = 0; i < 3; i++)
      {
         if(Dependencies[iE].size() <= i)
            continue;
         
         double Value = GetValue(Dependencies[iE][i]);
         if(Value < DependencyRanges[iE][i*2])
            Value = DependencyRanges[iE][i*2];
         if(Value > DependencyRanges[iE][i*2+1])
            Value = DependencyRanges[iE][i*2+1];
         V[i] = Value;
      }

      TF1 *Function = nullptr;
      
      if(Functions[iE] == nullptr)
      {
         if(Dependencies[iE].size() == 1)
            Function = new TF1(Form("Function%d", iE), (Formulas[iE] + "+0*x").c_str());
         if(Dependencies[iE].size() == 2)
            Function = new TF2(Form("Function%d", iE), (Formulas[iE] + "+0*x+0*y").c_str());
         if(Dependencies[iE].size() == 3)
            Function = new TF3(Form("Function%d", iE), (Formulas[iE] + "+0*x+0*y+0*z").c_str());
         if(Dependencies[iE].size() == 4)
            Function = new TF3(Form("Function%d", iE), (Formulas[iE] + "+0*x+0*y+0*z").c_str());

         Functions[iE] = Function;
      }
      else
         Function = Functions[iE];

      for(int i = 0; i < (int)Parameters[iE].size(); i++)
         Function->SetParameter(i, Parameters[iE][i]);
      if(Dependencies[iE].size() == 4)
         Function->SetParameter(Parameters[iE].size(), GetValue(Dependencies[iE][3]));
      double Result = Function->EvalPar(V);

      // cout << Formulas[iE] << endl;
      // cout << "P" << endl;
      // for(int i = 0; i < (int)Parameters[iE].size(); i++)
      //    cout << " " << Parameters[iE][i] << endl;
      // cout << "V" << endl;
      // cout << " " << V[0] << endl;
      // cout << " " << V[1] << endl;
      // cout << " " << V[2] << endl;
      // cout << Dependencies[iE].size() << endl;
      // cout << Function->EvalPar(V) << endl;

      return Result;
   }

   return -1;
}

double SingleJetCorrector::GetCorrectedPT()
{
   double Correction = GetCorrection();

   if(Correction < 0)
      return -1;

   return JetPT * Correction;
}

double SingleJetCorrector::GetValue(Type T)
{
   if(T == TypeNone)      return 0;
   if(T == TypeJetPT)     return JetPT;
   if(T == TypeJetEta)    return JetEta;
   if(T == TypeJetPhi)    return JetPhi;
   if(T == TypeJetArea)   return JetArea;
   if(T == TypeRho)       return Rho;

   return -1;
}

std::string SingleJetCorrector::Hack4(std::string Formula, char V, int N)
{
   int Size = Formula.size();
   for(int i = 0; i < Size; i++)
   {
      if(Formula[i] != V)
         continue;

      if(i != 0 && Formula[i-1] >= 'a' && Formula[i-1] <= 'z')   continue;
      if(i != 0 && Formula[i-1] >= 'A' && Formula[i-1] <= 'Z')   continue;
      if(i != Size - 1 && Formula[i+1] >= 'a' && Formula[i+1] <= 'z')   continue;
      if(i != Size - 1 && Formula[i+1] >= 'A' && Formula[i+1] <= 'Z')   continue;

      Formula.erase(i, 1);
      Formula.insert(i, Form("[%d]", N));
   }

   return Formula;
}

#endif
