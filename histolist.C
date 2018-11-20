void loopdir() {
  TFile *f1 = TFile::Open("output/output.root");
  TIter keyList(f1->GetListOfKeys());
  TKey *key;
  TCanvas c1; // CANVAS INSTANCE

  c1.Print("hsimple.ps[");
  while ((key = (TKey*)keyList())) {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      TH1 *h = (TH1*)key->ReadObj();
      h->Draw();
      c1.Print("hsimple.ps");
   }
  c1.Print("hsimple.ps]");
}













