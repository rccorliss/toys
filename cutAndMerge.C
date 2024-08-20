
void cutAndMerge(const char* directory, const char* outputFileName, const char* cut) {
    // Create a chain to hold the trees from multiple files
    TChain chain("residualtree");

    // Get the list of files in the directory
    TSystemDirectory dir(directory, directory);
    TList* files = dir.GetListOfFiles();

    // Loop over the files in the directory
    if (files) {
        TSystemFile* file;
        TString fileName;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fileName = file->GetName();
            if (!file->IsDirectory() && fileName.EndsWith(".root")) {
                // Open the file
                TFile* inputFile = TFile::Open(Form("%s/%s", directory, fileName.Data()), "READ");
                if (inputFile) {
                    // Get the tree from the file
                    TTree* tree = (TTree*)inputFile->Get("residualtree");
                    if (tree) {
                        // Apply the cut
                        TCut cutObj(cut);
                        tree->Draw(">>entryList", cutObj);
                        TEntryList* entryList = (TEntryList*)gDirectory->Get("entryList");
                        tree->SetEntryList(entryList);

                        // Add the tree to the chain
                        chain.Add(Form("%s/%s", directory, fileName.Data()));
                    }
                    inputFile->Close();
                }
            }
        }
    }

    // Merge the trees in the chain
    TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
    TTree* mergedTree = chain.MergeTrees();
    mergedTree->SetName("residualtree");
    mergedTree->Write();
    outputFile->Close();
}
