// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME RootDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/interface/Packager.h"
#include "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/interface/WSTFileWrapper.h"
#include "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/interface/LinearInterp.h"
#include "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/interface/InitialFit.h"
#include "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/interface/RooGaussBern2D.h"
#include "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/interface/ReplacementMap.h"
#include "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/interface/Normalization_8TeV.h"
#include "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/interface/Normalization_13TeV.h"
#include "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/interface/SimultaneousFit.h"
#include "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/interface/LCRooChi2Var.h"
#include "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/interface/FinalModelConstruction.h"
#include "/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/interface/LCRooAddition.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *Normalization_8TeV_Dictionary();
   static void Normalization_8TeV_TClassManip(TClass*);
   static void *new_Normalization_8TeV(void *p = 0);
   static void *newArray_Normalization_8TeV(Long_t size, void *p);
   static void delete_Normalization_8TeV(void *p);
   static void deleteArray_Normalization_8TeV(void *p);
   static void destruct_Normalization_8TeV(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Normalization_8TeV*)
   {
      ::Normalization_8TeV *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Normalization_8TeV));
      static ::ROOT::TGenericClassInfo 
         instance("Normalization_8TeV", "interface/Normalization_8TeV.h", 19,
                  typeid(::Normalization_8TeV), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Normalization_8TeV_Dictionary, isa_proxy, 0,
                  sizeof(::Normalization_8TeV) );
      instance.SetNew(&new_Normalization_8TeV);
      instance.SetNewArray(&newArray_Normalization_8TeV);
      instance.SetDelete(&delete_Normalization_8TeV);
      instance.SetDeleteArray(&deleteArray_Normalization_8TeV);
      instance.SetDestructor(&destruct_Normalization_8TeV);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Normalization_8TeV*)
   {
      return GenerateInitInstanceLocal((::Normalization_8TeV*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Normalization_8TeV*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Normalization_8TeV_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Normalization_8TeV*)0x0)->GetClass();
      Normalization_8TeV_TClassManip(theClass);
   return theClass;
   }

   static void Normalization_8TeV_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *WSTFileWrapper_Dictionary();
   static void WSTFileWrapper_TClassManip(TClass*);
   static void delete_WSTFileWrapper(void *p);
   static void deleteArray_WSTFileWrapper(void *p);
   static void destruct_WSTFileWrapper(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WSTFileWrapper*)
   {
      ::WSTFileWrapper *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::WSTFileWrapper));
      static ::ROOT::TGenericClassInfo 
         instance("WSTFileWrapper", "interface/WSTFileWrapper.h", 12,
                  typeid(::WSTFileWrapper), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &WSTFileWrapper_Dictionary, isa_proxy, 0,
                  sizeof(::WSTFileWrapper) );
      instance.SetDelete(&delete_WSTFileWrapper);
      instance.SetDeleteArray(&deleteArray_WSTFileWrapper);
      instance.SetDestructor(&destruct_WSTFileWrapper);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WSTFileWrapper*)
   {
      return GenerateInitInstanceLocal((::WSTFileWrapper*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::WSTFileWrapper*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *WSTFileWrapper_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::WSTFileWrapper*)0x0)->GetClass();
      WSTFileWrapper_TClassManip(theClass);
   return theClass;
   }

   static void WSTFileWrapper_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Packager_Dictionary();
   static void Packager_TClassManip(TClass*);
   static void delete_Packager(void *p);
   static void deleteArray_Packager(void *p);
   static void destruct_Packager(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Packager*)
   {
      ::Packager *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Packager));
      static ::ROOT::TGenericClassInfo 
         instance("Packager", "interface/Packager.h", 17,
                  typeid(::Packager), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Packager_Dictionary, isa_proxy, 0,
                  sizeof(::Packager) );
      instance.SetDelete(&delete_Packager);
      instance.SetDeleteArray(&deleteArray_Packager);
      instance.SetDestructor(&destruct_Packager);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Packager*)
   {
      return GenerateInitInstanceLocal((::Packager*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Packager*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Packager_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Packager*)0x0)->GetClass();
      Packager_TClassManip(theClass);
   return theClass;
   }

   static void Packager_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *LinearInterp_Dictionary();
   static void LinearInterp_TClassManip(TClass*);
   static void delete_LinearInterp(void *p);
   static void deleteArray_LinearInterp(void *p);
   static void destruct_LinearInterp(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LinearInterp*)
   {
      ::LinearInterp *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::LinearInterp));
      static ::ROOT::TGenericClassInfo 
         instance("LinearInterp", "interface/LinearInterp.h", 15,
                  typeid(::LinearInterp), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &LinearInterp_Dictionary, isa_proxy, 0,
                  sizeof(::LinearInterp) );
      instance.SetDelete(&delete_LinearInterp);
      instance.SetDeleteArray(&deleteArray_LinearInterp);
      instance.SetDestructor(&destruct_LinearInterp);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LinearInterp*)
   {
      return GenerateInitInstanceLocal((::LinearInterp*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LinearInterp*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *LinearInterp_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::LinearInterp*)0x0)->GetClass();
      LinearInterp_TClassManip(theClass);
   return theClass;
   }

   static void LinearInterp_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *InitialFit_Dictionary();
   static void InitialFit_TClassManip(TClass*);
   static void delete_InitialFit(void *p);
   static void deleteArray_InitialFit(void *p);
   static void destruct_InitialFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::InitialFit*)
   {
      ::InitialFit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::InitialFit));
      static ::ROOT::TGenericClassInfo 
         instance("InitialFit", "interface/InitialFit.h", 17,
                  typeid(::InitialFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &InitialFit_Dictionary, isa_proxy, 0,
                  sizeof(::InitialFit) );
      instance.SetDelete(&delete_InitialFit);
      instance.SetDeleteArray(&deleteArray_InitialFit);
      instance.SetDestructor(&destruct_InitialFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::InitialFit*)
   {
      return GenerateInitInstanceLocal((::InitialFit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::InitialFit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *InitialFit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::InitialFit*)0x0)->GetClass();
      InitialFit_TClassManip(theClass);
   return theClass;
   }

   static void InitialFit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_RooGaussBern2D(void *p = 0);
   static void *newArray_RooGaussBern2D(Long_t size, void *p);
   static void delete_RooGaussBern2D(void *p);
   static void deleteArray_RooGaussBern2D(void *p);
   static void destruct_RooGaussBern2D(void *p);
   static void streamer_RooGaussBern2D(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooGaussBern2D*)
   {
      ::RooGaussBern2D *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooGaussBern2D >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooGaussBern2D", ::RooGaussBern2D::Class_Version(), "interface/RooGaussBern2D.h", 12,
                  typeid(::RooGaussBern2D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooGaussBern2D::Dictionary, isa_proxy, 16,
                  sizeof(::RooGaussBern2D) );
      instance.SetNew(&new_RooGaussBern2D);
      instance.SetNewArray(&newArray_RooGaussBern2D);
      instance.SetDelete(&delete_RooGaussBern2D);
      instance.SetDeleteArray(&deleteArray_RooGaussBern2D);
      instance.SetDestructor(&destruct_RooGaussBern2D);
      instance.SetStreamerFunc(&streamer_RooGaussBern2D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooGaussBern2D*)
   {
      return GenerateInitInstanceLocal((::RooGaussBern2D*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooGaussBern2D*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *ReplacementMap_Dictionary();
   static void ReplacementMap_TClassManip(TClass*);
   static void delete_ReplacementMap(void *p);
   static void deleteArray_ReplacementMap(void *p);
   static void destruct_ReplacementMap(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ReplacementMap*)
   {
      ::ReplacementMap *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ReplacementMap));
      static ::ROOT::TGenericClassInfo 
         instance("ReplacementMap", "interface/ReplacementMap.h", 16,
                  typeid(::ReplacementMap), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ReplacementMap_Dictionary, isa_proxy, 0,
                  sizeof(::ReplacementMap) );
      instance.SetDelete(&delete_ReplacementMap);
      instance.SetDeleteArray(&deleteArray_ReplacementMap);
      instance.SetDestructor(&destruct_ReplacementMap);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ReplacementMap*)
   {
      return GenerateInitInstanceLocal((::ReplacementMap*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ReplacementMap*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ReplacementMap_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ReplacementMap*)0x0)->GetClass();
      ReplacementMap_TClassManip(theClass);
   return theClass;
   }

   static void ReplacementMap_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Normalization_13TeV_Dictionary();
   static void Normalization_13TeV_TClassManip(TClass*);
   static void *new_Normalization_13TeV(void *p = 0);
   static void *newArray_Normalization_13TeV(Long_t size, void *p);
   static void delete_Normalization_13TeV(void *p);
   static void deleteArray_Normalization_13TeV(void *p);
   static void destruct_Normalization_13TeV(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Normalization_13TeV*)
   {
      ::Normalization_13TeV *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Normalization_13TeV));
      static ::ROOT::TGenericClassInfo 
         instance("Normalization_13TeV", "interface/Normalization_13TeV.h", 19,
                  typeid(::Normalization_13TeV), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Normalization_13TeV_Dictionary, isa_proxy, 0,
                  sizeof(::Normalization_13TeV) );
      instance.SetNew(&new_Normalization_13TeV);
      instance.SetNewArray(&newArray_Normalization_13TeV);
      instance.SetDelete(&delete_Normalization_13TeV);
      instance.SetDeleteArray(&deleteArray_Normalization_13TeV);
      instance.SetDestructor(&destruct_Normalization_13TeV);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Normalization_13TeV*)
   {
      return GenerateInitInstanceLocal((::Normalization_13TeV*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Normalization_13TeV*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Normalization_13TeV_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::Normalization_13TeV*)0x0)->GetClass();
      Normalization_13TeV_TClassManip(theClass);
   return theClass;
   }

   static void Normalization_13TeV_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *SimultaneousFit_Dictionary();
   static void SimultaneousFit_TClassManip(TClass*);
   static void delete_SimultaneousFit(void *p);
   static void deleteArray_SimultaneousFit(void *p);
   static void destruct_SimultaneousFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SimultaneousFit*)
   {
      ::SimultaneousFit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SimultaneousFit));
      static ::ROOT::TGenericClassInfo 
         instance("SimultaneousFit", "interface/SimultaneousFit.h", 20,
                  typeid(::SimultaneousFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SimultaneousFit_Dictionary, isa_proxy, 0,
                  sizeof(::SimultaneousFit) );
      instance.SetDelete(&delete_SimultaneousFit);
      instance.SetDeleteArray(&deleteArray_SimultaneousFit);
      instance.SetDestructor(&destruct_SimultaneousFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SimultaneousFit*)
   {
      return GenerateInitInstanceLocal((::SimultaneousFit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SimultaneousFit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SimultaneousFit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SimultaneousFit*)0x0)->GetClass();
      SimultaneousFit_TClassManip(theClass);
   return theClass;
   }

   static void SimultaneousFit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void delete_LCRooChi2Var(void *p);
   static void deleteArray_LCRooChi2Var(void *p);
   static void destruct_LCRooChi2Var(void *p);
   static void streamer_LCRooChi2Var(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LCRooChi2Var*)
   {
      ::LCRooChi2Var *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LCRooChi2Var >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LCRooChi2Var", ::LCRooChi2Var::Class_Version(), "interface/LCRooChi2Var.h", 24,
                  typeid(::LCRooChi2Var), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LCRooChi2Var::Dictionary, isa_proxy, 16,
                  sizeof(::LCRooChi2Var) );
      instance.SetDelete(&delete_LCRooChi2Var);
      instance.SetDeleteArray(&deleteArray_LCRooChi2Var);
      instance.SetDestructor(&destruct_LCRooChi2Var);
      instance.SetStreamerFunc(&streamer_LCRooChi2Var);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LCRooChi2Var*)
   {
      return GenerateInitInstanceLocal((::LCRooChi2Var*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LCRooChi2Var*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *FinalModelConstruction_Dictionary();
   static void FinalModelConstruction_TClassManip(TClass*);
   static void delete_FinalModelConstruction(void *p);
   static void deleteArray_FinalModelConstruction(void *p);
   static void destruct_FinalModelConstruction(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FinalModelConstruction*)
   {
      ::FinalModelConstruction *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::FinalModelConstruction));
      static ::ROOT::TGenericClassInfo 
         instance("FinalModelConstruction", "interface/FinalModelConstruction.h", 26,
                  typeid(::FinalModelConstruction), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &FinalModelConstruction_Dictionary, isa_proxy, 0,
                  sizeof(::FinalModelConstruction) );
      instance.SetDelete(&delete_FinalModelConstruction);
      instance.SetDeleteArray(&deleteArray_FinalModelConstruction);
      instance.SetDestructor(&destruct_FinalModelConstruction);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FinalModelConstruction*)
   {
      return GenerateInitInstanceLocal((::FinalModelConstruction*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::FinalModelConstruction*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *FinalModelConstruction_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::FinalModelConstruction*)0x0)->GetClass();
      FinalModelConstruction_TClassManip(theClass);
   return theClass;
   }

   static void FinalModelConstruction_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_LCRooAddition(void *p = 0);
   static void *newArray_LCRooAddition(Long_t size, void *p);
   static void delete_LCRooAddition(void *p);
   static void deleteArray_LCRooAddition(void *p);
   static void destruct_LCRooAddition(void *p);
   static void streamer_LCRooAddition(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LCRooAddition*)
   {
      ::LCRooAddition *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LCRooAddition >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LCRooAddition", ::LCRooAddition::Class_Version(), "interface/LCRooAddition.h", 29,
                  typeid(::LCRooAddition), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LCRooAddition::Dictionary, isa_proxy, 16,
                  sizeof(::LCRooAddition) );
      instance.SetNew(&new_LCRooAddition);
      instance.SetNewArray(&newArray_LCRooAddition);
      instance.SetDelete(&delete_LCRooAddition);
      instance.SetDeleteArray(&deleteArray_LCRooAddition);
      instance.SetDestructor(&destruct_LCRooAddition);
      instance.SetStreamerFunc(&streamer_LCRooAddition);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LCRooAddition*)
   {
      return GenerateInitInstanceLocal((::LCRooAddition*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::LCRooAddition*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RooGaussBern2D::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooGaussBern2D::Class_Name()
{
   return "RooGaussBern2D";
}

//______________________________________________________________________________
const char *RooGaussBern2D::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooGaussBern2D*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooGaussBern2D::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooGaussBern2D*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooGaussBern2D::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooGaussBern2D*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooGaussBern2D::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooGaussBern2D*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr LCRooChi2Var::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LCRooChi2Var::Class_Name()
{
   return "LCRooChi2Var";
}

//______________________________________________________________________________
const char *LCRooChi2Var::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LCRooChi2Var*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LCRooChi2Var::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LCRooChi2Var*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LCRooChi2Var::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LCRooChi2Var*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LCRooChi2Var::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LCRooChi2Var*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr LCRooAddition::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LCRooAddition::Class_Name()
{
   return "LCRooAddition";
}

//______________________________________________________________________________
const char *LCRooAddition::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LCRooAddition*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LCRooAddition::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LCRooAddition*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LCRooAddition::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LCRooAddition*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LCRooAddition::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LCRooAddition*)0x0)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Normalization_8TeV(void *p) {
      return  p ? new(p) ::Normalization_8TeV : new ::Normalization_8TeV;
   }
   static void *newArray_Normalization_8TeV(Long_t nElements, void *p) {
      return p ? new(p) ::Normalization_8TeV[nElements] : new ::Normalization_8TeV[nElements];
   }
   // Wrapper around operator delete
   static void delete_Normalization_8TeV(void *p) {
      delete ((::Normalization_8TeV*)p);
   }
   static void deleteArray_Normalization_8TeV(void *p) {
      delete [] ((::Normalization_8TeV*)p);
   }
   static void destruct_Normalization_8TeV(void *p) {
      typedef ::Normalization_8TeV current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Normalization_8TeV

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WSTFileWrapper(void *p) {
      delete ((::WSTFileWrapper*)p);
   }
   static void deleteArray_WSTFileWrapper(void *p) {
      delete [] ((::WSTFileWrapper*)p);
   }
   static void destruct_WSTFileWrapper(void *p) {
      typedef ::WSTFileWrapper current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WSTFileWrapper

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Packager(void *p) {
      delete ((::Packager*)p);
   }
   static void deleteArray_Packager(void *p) {
      delete [] ((::Packager*)p);
   }
   static void destruct_Packager(void *p) {
      typedef ::Packager current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Packager

namespace ROOT {
   // Wrapper around operator delete
   static void delete_LinearInterp(void *p) {
      delete ((::LinearInterp*)p);
   }
   static void deleteArray_LinearInterp(void *p) {
      delete [] ((::LinearInterp*)p);
   }
   static void destruct_LinearInterp(void *p) {
      typedef ::LinearInterp current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LinearInterp

namespace ROOT {
   // Wrapper around operator delete
   static void delete_InitialFit(void *p) {
      delete ((::InitialFit*)p);
   }
   static void deleteArray_InitialFit(void *p) {
      delete [] ((::InitialFit*)p);
   }
   static void destruct_InitialFit(void *p) {
      typedef ::InitialFit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::InitialFit

//______________________________________________________________________________
void RooGaussBern2D::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooGaussBern2D.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      x.Streamer(R__b);
      y.Streamer(R__b);
      polParamsMean.Streamer(R__b);
      polParamsSigma.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, RooGaussBern2D::IsA());
   } else {
      R__c = R__b.WriteVersion(RooGaussBern2D::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      x.Streamer(R__b);
      y.Streamer(R__b);
      polParamsMean.Streamer(R__b);
      polParamsSigma.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RooGaussBern2D(void *p) {
      return  p ? new(p) ::RooGaussBern2D : new ::RooGaussBern2D;
   }
   static void *newArray_RooGaussBern2D(Long_t nElements, void *p) {
      return p ? new(p) ::RooGaussBern2D[nElements] : new ::RooGaussBern2D[nElements];
   }
   // Wrapper around operator delete
   static void delete_RooGaussBern2D(void *p) {
      delete ((::RooGaussBern2D*)p);
   }
   static void deleteArray_RooGaussBern2D(void *p) {
      delete [] ((::RooGaussBern2D*)p);
   }
   static void destruct_RooGaussBern2D(void *p) {
      typedef ::RooGaussBern2D current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooGaussBern2D(TBuffer &buf, void *obj) {
      ((::RooGaussBern2D*)obj)->::RooGaussBern2D::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooGaussBern2D

namespace ROOT {
   // Wrapper around operator delete
   static void delete_ReplacementMap(void *p) {
      delete ((::ReplacementMap*)p);
   }
   static void deleteArray_ReplacementMap(void *p) {
      delete [] ((::ReplacementMap*)p);
   }
   static void destruct_ReplacementMap(void *p) {
      typedef ::ReplacementMap current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ReplacementMap

namespace ROOT {
   // Wrappers around operator new
   static void *new_Normalization_13TeV(void *p) {
      return  p ? new(p) ::Normalization_13TeV : new ::Normalization_13TeV;
   }
   static void *newArray_Normalization_13TeV(Long_t nElements, void *p) {
      return p ? new(p) ::Normalization_13TeV[nElements] : new ::Normalization_13TeV[nElements];
   }
   // Wrapper around operator delete
   static void delete_Normalization_13TeV(void *p) {
      delete ((::Normalization_13TeV*)p);
   }
   static void deleteArray_Normalization_13TeV(void *p) {
      delete [] ((::Normalization_13TeV*)p);
   }
   static void destruct_Normalization_13TeV(void *p) {
      typedef ::Normalization_13TeV current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Normalization_13TeV

namespace ROOT {
   // Wrapper around operator delete
   static void delete_SimultaneousFit(void *p) {
      delete ((::SimultaneousFit*)p);
   }
   static void deleteArray_SimultaneousFit(void *p) {
      delete [] ((::SimultaneousFit*)p);
   }
   static void destruct_SimultaneousFit(void *p) {
      typedef ::SimultaneousFit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SimultaneousFit

//______________________________________________________________________________
void LCRooChi2Var::Streamer(TBuffer &R__b)
{
   // Stream an object of class LCRooChi2Var.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsOptTestStatistic::Streamer(R__b);
      R__b >> _data;
      void *ptr__etype = (void*)&_etype;
      R__b >> *reinterpret_cast<Int_t*>(ptr__etype);
      void *ptr__funcMode = (void*)&_funcMode;
      R__b >> *reinterpret_cast<Int_t*>(ptr__funcMode);
      R__b.CheckByteCount(R__s, R__c, LCRooChi2Var::IsA());
   } else {
      R__c = R__b.WriteVersion(LCRooChi2Var::IsA(), kTRUE);
      RooAbsOptTestStatistic::Streamer(R__b);
      R__b << _data;
      R__b << (Int_t)_etype;
      R__b << (Int_t)_funcMode;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_LCRooChi2Var(void *p) {
      delete ((::LCRooChi2Var*)p);
   }
   static void deleteArray_LCRooChi2Var(void *p) {
      delete [] ((::LCRooChi2Var*)p);
   }
   static void destruct_LCRooChi2Var(void *p) {
      typedef ::LCRooChi2Var current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_LCRooChi2Var(TBuffer &buf, void *obj) {
      ((::LCRooChi2Var*)obj)->::LCRooChi2Var::Streamer(buf);
   }
} // end of namespace ROOT for class ::LCRooChi2Var

namespace ROOT {
   // Wrapper around operator delete
   static void delete_FinalModelConstruction(void *p) {
      delete ((::FinalModelConstruction*)p);
   }
   static void deleteArray_FinalModelConstruction(void *p) {
      delete [] ((::FinalModelConstruction*)p);
   }
   static void destruct_FinalModelConstruction(void *p) {
      typedef ::FinalModelConstruction current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::FinalModelConstruction

//______________________________________________________________________________
void LCRooAddition::Streamer(TBuffer &R__b)
{
   // Stream an object of class LCRooAddition.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsReal::Streamer(R__b);
      R__b >> _pdf;
      {
         map<int,RooDataHist*> &R__stl =  _datasets;
         R__stl.clear();
         TClass *R__tcl2 = TBuffer::GetClass(typeid(class RooDataHist *));
         if (R__tcl2==0) {
            Error("_datasets streamer","Missing the TClass object for class RooDataHist *!");
            return;
         }
         int R__i, R__n;
         R__b >> R__n;
         for (R__i = 0; R__i < R__n; R__i++) {
            int R__t;
            R__b >> R__t;
            RooDataHist* R__t2;
            R__t2 = (RooDataHist*)R__b.ReadObjectAny(R__tcl2);
            typedef int Value_t;
            std::pair<Value_t const, class RooDataHist * > R__t3(R__t,R__t2);
            R__stl.insert(R__t3);
         }
      }
      R__b >> _MH;
      R__b >> _mgg;
      {
         map<int,LCRooChi2Var*> &R__stl =  _chi2map;
         R__stl.clear();
         TClass *R__tcl2 = TBuffer::GetClass(typeid(class LCRooChi2Var *));
         if (R__tcl2==0) {
            Error("_chi2map streamer","Missing the TClass object for class LCRooChi2Var *!");
            return;
         }
         int R__i, R__n;
         R__b >> R__n;
         for (R__i = 0; R__i < R__n; R__i++) {
            int R__t;
            R__b >> R__t;
            LCRooChi2Var* R__t2;
            R__t2 = (LCRooChi2Var*)R__b.ReadObjectAny(R__tcl2);
            typedef int Value_t;
            std::pair<Value_t const, class LCRooChi2Var * > R__t3(R__t,R__t2);
            R__stl.insert(R__t3);
         }
      }
      _ownedList.Streamer(R__b);
      _set.Streamer(R__b);
      _cacheMgr.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, LCRooAddition::IsA());
   } else {
      R__c = R__b.WriteVersion(LCRooAddition::IsA(), kTRUE);
      RooAbsReal::Streamer(R__b);
      R__b << _pdf;
      {
         map<int,RooDataHist*> &R__stl =  _datasets;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            map<int,RooDataHist*>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << ((*R__k).first );
            R__b << ((*R__k).second);
            }
         }
      }
      R__b << _MH;
      R__b << _mgg;
      {
         map<int,LCRooChi2Var*> &R__stl =  _chi2map;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            map<int,LCRooChi2Var*>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << ((*R__k).first );
            R__b << ((*R__k).second);
            }
         }
      }
      _ownedList.Streamer(R__b);
      _set.Streamer(R__b);
      _cacheMgr.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LCRooAddition(void *p) {
      return  p ? new(p) ::LCRooAddition : new ::LCRooAddition;
   }
   static void *newArray_LCRooAddition(Long_t nElements, void *p) {
      return p ? new(p) ::LCRooAddition[nElements] : new ::LCRooAddition[nElements];
   }
   // Wrapper around operator delete
   static void delete_LCRooAddition(void *p) {
      delete ((::LCRooAddition*)p);
   }
   static void deleteArray_LCRooAddition(void *p) {
      delete [] ((::LCRooAddition*)p);
   }
   static void destruct_LCRooAddition(void *p) {
      typedef ::LCRooAddition current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_LCRooAddition(TBuffer &buf, void *obj) {
      ((::LCRooAddition*)obj)->::LCRooAddition::Streamer(buf);
   }
} // end of namespace ROOT for class ::LCRooAddition

namespace ROOT {
   static TClass *maplEintcORooDataHistmUgR_Dictionary();
   static void maplEintcORooDataHistmUgR_TClassManip(TClass*);
   static void *new_maplEintcORooDataHistmUgR(void *p = 0);
   static void *newArray_maplEintcORooDataHistmUgR(Long_t size, void *p);
   static void delete_maplEintcORooDataHistmUgR(void *p);
   static void deleteArray_maplEintcORooDataHistmUgR(void *p);
   static void destruct_maplEintcORooDataHistmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,RooDataHist*>*)
   {
      map<int,RooDataHist*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,RooDataHist*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,RooDataHist*>", -2, "map", 99,
                  typeid(map<int,RooDataHist*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcORooDataHistmUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,RooDataHist*>) );
      instance.SetNew(&new_maplEintcORooDataHistmUgR);
      instance.SetNewArray(&newArray_maplEintcORooDataHistmUgR);
      instance.SetDelete(&delete_maplEintcORooDataHistmUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcORooDataHistmUgR);
      instance.SetDestructor(&destruct_maplEintcORooDataHistmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,RooDataHist*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,RooDataHist*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcORooDataHistmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,RooDataHist*>*)0x0)->GetClass();
      maplEintcORooDataHistmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcORooDataHistmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcORooDataHistmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,RooDataHist*> : new map<int,RooDataHist*>;
   }
   static void *newArray_maplEintcORooDataHistmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,RooDataHist*>[nElements] : new map<int,RooDataHist*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcORooDataHistmUgR(void *p) {
      delete ((map<int,RooDataHist*>*)p);
   }
   static void deleteArray_maplEintcORooDataHistmUgR(void *p) {
      delete [] ((map<int,RooDataHist*>*)p);
   }
   static void destruct_maplEintcORooDataHistmUgR(void *p) {
      typedef map<int,RooDataHist*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,RooDataHist*>

namespace ROOT {
   static TClass *maplEintcOLCRooChi2VarmUgR_Dictionary();
   static void maplEintcOLCRooChi2VarmUgR_TClassManip(TClass*);
   static void *new_maplEintcOLCRooChi2VarmUgR(void *p = 0);
   static void *newArray_maplEintcOLCRooChi2VarmUgR(Long_t size, void *p);
   static void delete_maplEintcOLCRooChi2VarmUgR(void *p);
   static void deleteArray_maplEintcOLCRooChi2VarmUgR(void *p);
   static void destruct_maplEintcOLCRooChi2VarmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,LCRooChi2Var*>*)
   {
      map<int,LCRooChi2Var*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,LCRooChi2Var*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,LCRooChi2Var*>", -2, "map", 99,
                  typeid(map<int,LCRooChi2Var*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOLCRooChi2VarmUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,LCRooChi2Var*>) );
      instance.SetNew(&new_maplEintcOLCRooChi2VarmUgR);
      instance.SetNewArray(&newArray_maplEintcOLCRooChi2VarmUgR);
      instance.SetDelete(&delete_maplEintcOLCRooChi2VarmUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOLCRooChi2VarmUgR);
      instance.SetDestructor(&destruct_maplEintcOLCRooChi2VarmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,LCRooChi2Var*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,LCRooChi2Var*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOLCRooChi2VarmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,LCRooChi2Var*>*)0x0)->GetClass();
      maplEintcOLCRooChi2VarmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOLCRooChi2VarmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOLCRooChi2VarmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,LCRooChi2Var*> : new map<int,LCRooChi2Var*>;
   }
   static void *newArray_maplEintcOLCRooChi2VarmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,LCRooChi2Var*>[nElements] : new map<int,LCRooChi2Var*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOLCRooChi2VarmUgR(void *p) {
      delete ((map<int,LCRooChi2Var*>*)p);
   }
   static void deleteArray_maplEintcOLCRooChi2VarmUgR(void *p) {
      delete [] ((map<int,LCRooChi2Var*>*)p);
   }
   static void destruct_maplEintcOLCRooChi2VarmUgR(void *p) {
      typedef map<int,LCRooChi2Var*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,LCRooChi2Var*>

namespace {
  void TriggerDictionaryInitialization_RootDict_Impl() {
    static const char* headers[] = {
"interface/Packager.h",
"interface/WSTFileWrapper.h",
"interface/LinearInterp.h",
"interface/InitialFit.h",
"interface/RooGaussBern2D.h",
"interface/ReplacementMap.h",
"interface/Normalization_8TeV.h",
"interface/Normalization_13TeV.h",
"interface/SimultaneousFit.h",
"interface/LCRooChi2Var.h",
"interface/FinalModelConstruction.h",
"interface/LCRooAddition.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.12.07-gnimlf5//include",
"/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src",
"/cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.12.07-gnimlf5/include",
"/afs/cern.ch/work/z/zewang/private/HZGamma/flashggfinalfit/CMSSW_10_2_13/src/flashggFinalFit/Signal/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "RootDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$interface/Packager.h")))  Normalization_8TeV;
class __attribute__((annotate("$clingAutoload$interface/Packager.h")))  WSTFileWrapper;
class __attribute__((annotate("$clingAutoload$interface/Packager.h")))  Packager;
class __attribute__((annotate("$clingAutoload$interface/LinearInterp.h")))  LinearInterp;
class __attribute__((annotate("$clingAutoload$interface/InitialFit.h")))  InitialFit;
class __attribute__((annotate("$clingAutoload$interface/RooGaussBern2D.h")))  RooGaussBern2D;
class __attribute__((annotate("$clingAutoload$interface/ReplacementMap.h")))  ReplacementMap;
class __attribute__((annotate("$clingAutoload$interface/Normalization_13TeV.h")))  Normalization_13TeV;
class __attribute__((annotate("$clingAutoload$interface/SimultaneousFit.h")))  SimultaneousFit;
class __attribute__((annotate(R"ATTRDUMP(Chi^2 function of p.d.f w.r.t a binned dataset)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/LCRooChi2Var.h")))  LCRooChi2Var;
class __attribute__((annotate("$clingAutoload$interface/FinalModelConstruction.h")))  FinalModelConstruction;
class __attribute__((annotate(R"ATTRDUMP(Sum of RooAbsReal objects)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$interface/LCRooAddition.h")))  LCRooAddition;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "RootDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "interface/Packager.h"
#include "interface/WSTFileWrapper.h"
#include "interface/LinearInterp.h"
#include "interface/InitialFit.h"
#include "interface/RooGaussBern2D.h"
#include "interface/ReplacementMap.h"
#include "interface/Normalization_8TeV.h"
#include "interface/Normalization_13TeV.h"
#include "interface/SimultaneousFit.h"
#include "interface/LCRooChi2Var.h"
#include "interface/FinalModelConstruction.h"
#include "interface/LCRooAddition.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"FinalModelConstruction", payloadCode, "@",
"InitialFit", payloadCode, "@",
"LCRooAddition", payloadCode, "@",
"LCRooChi2Var", payloadCode, "@",
"LinearInterp", payloadCode, "@",
"Normalization_13TeV", payloadCode, "@",
"Normalization_8TeV", payloadCode, "@",
"Packager", payloadCode, "@",
"ReplacementMap", payloadCode, "@",
"RooGaussBern2D", payloadCode, "@",
"SimultaneousFit", payloadCode, "@",
"WSTFileWrapper", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RootDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RootDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RootDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RootDict() {
  TriggerDictionaryInitialization_RootDict_Impl();
}
