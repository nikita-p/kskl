{
    gROOT->ProcessLine(".L handler/fit.cpp+");
    gROOT->ProcessLine(".L handler/dataset.cpp+");

    gROOT->ProcessLine(" Dataset d(\"Outputs/trees/model\", \"Outputs/trees/2011\", \"Inputs/2011/lum.dat\", \"Inputs/2011/rad_cors.dat\")");
}
