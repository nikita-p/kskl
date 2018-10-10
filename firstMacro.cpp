{
    gROOT->ProcessLine(".L readingTree/events.cpp+");
    gROOT->ProcessLine(".L readingTree/Trph.C+");
    cout << "Year: ";
    string year;
    cin >> year;
    events(year);
}
