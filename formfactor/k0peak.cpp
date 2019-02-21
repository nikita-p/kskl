{
//========= Macro generated from object: Graph/Graph
//========= by ROOT version6.14/04
   
   Double_t Graph_fx1001[25] = {
   1.00407,
   1.01047,
   1.01296,
   1.01507,
   1.01611,
   1.01716,
   1.01716,
   1.01805,
   1.01912,
   1.01921,
   1.01942,
   1.0199,
   1.02122,
   1.02131,
   1.02208,
   1.02274,
   1.02326,
   1.02532,
   1.02796,
   1.02909,
   1.03391,
   1.04003,
   1.04986,
   1.05086,
   1.05995};
   Double_t Graph_fy1001[25] = {
   6.87,
   42.16,
   96.74,
   219.53,
   366.33,
   628.15,
   624.76,
   996.62,
   1413.65,
   1433.05,
   1434.84,
   1341.91,
   833.2,
   807.54,
   582.93,
   443.71,
   377.77,
   199.26,
   115.93,
   96.96,
   50.12,
   31.27,
   16.93,
   17.47,
   12.09};
   Double_t Graph_fex1001[25] = {
   8e-06,
   1e-05,
   7e-06,
   1.2e-05,
   1e-05,
   1.2e-05,
   1.3e-05,
   2.1e-05,
   1.6e-05,
   1.9e-05,
   2.8e-05,
   1.2e-05,
   2.1e-05,
   9e-06,
   2.1e-05,
   1.9e-05,
   2.5e-05,
   3.1e-05,
   1.5e-05,
   1.4e-05,
   1.1e-05,
   3.5e-05,
   1.1e-05,
   3.1e-05,
   1.5e-05};
   Double_t Graph_fey1001[25] = {
   0.42,
   0.47,
   1,
   5.02,
   3.33,
   2.95,
   9.89,
   4.28,
   6.02,
   15.03,
   18.4,
   4.74,
   4.89,
   10.36,
   4.03,
   4.38,
   5.31,
   4.97,
   1.7,
   3,
   1.26,
   1.01,
   0.5,
   0.94,
   0.71};
   TGraphErrors *gre = new TGraphErrors(25,Graph_fx1001,Graph_fy1001,Graph_fex1001,Graph_fey1001);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);
   gStyle->SetOptFit(11111);
   
   TH1F *Graph_Graph1001 = new TH1F("gre","Graph",100,0.998468,1.06555);
   Graph_Graph1001->SetMinimum(5.805);
   Graph_Graph1001->SetMaximum(1597.92);
   Graph_Graph1001->SetDirectory(0);
   Graph_Graph1001->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph1001->SetLineColor(ci);
   Graph_Graph1001->GetXaxis()->SetLabelFont(42);
   Graph_Graph1001->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1001->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1001->GetXaxis()->SetTitleFont(42);
   Graph_Graph1001->GetYaxis()->SetLabelFont(42);
   Graph_Graph1001->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1001->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1001->GetYaxis()->SetTitleOffset(0);
   Graph_Graph1001->GetYaxis()->SetTitleFont(42);
   Graph_Graph1001->GetZaxis()->SetLabelFont(42);
   Graph_Graph1001->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1001->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1001);
   
   gre->Draw("");
}
