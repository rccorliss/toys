bool checkAxisInterpolationBounds(float val, TAxis* axis){
    //get the middle point of the first bin on the axis:
    float lowerBound = axis->GetBinLowEdge(1)+axis->GetBinWidth(1)/2;

    //get the lower edge of the last bin
    int lastbin=axis->GetNbins();
    float upperBound = axis->GetBinLowEdge(lastbin)+axis->GetBinWidth(lastbin)/2;

    if (val<lowerBound || val>upperBound) {
        printf("value out of interpolation bounds: %f not within [%f,%f]\n",val,lowerBound,upperBound);
        return false;
    }

    return true;
}

bool checkBounds(float phi, float r, float z, TH3* h){
    if (!checkAxisInterpolationBounds(phi, h->GetXaxis())) {
        printf("failed phi bounds check\n");
        return false;
    }
    if (!checkAxisInterpolationBounds(r, h->GetYaxis())) {
        printf("failed r bounds check\n");
        return false;
    }
    if (!checkAxisInterpolationBounds(z, h->GetZaxis())) {
        printf("failed z bounds check\n");
        return false;
    }
    return true;
}

TVector3 correctPosition(TVector3 pos, TH3* hDPint, TH3* hDRint, TH3* hDZint, bool isRadians){
    //printf("histogram: %s pos: %f %f %f\n", pos.Phi(), pos.Perp(), pos.Z(),hDPint->GetName());
    //printf("second bin midpoints: %f %f %f\n", hDPint->GetXaxis()->GetBinCenter(2), hDPint->GetYaxis()->GetBinCenter(2), hDPint->GetZaxis()->GetBinCenter(2));
    //printf("second to last bin midpoints: %f %f %f\n", hDPint->GetXaxis()->GetBinCenter(hDPint->GetNbinsX()-1), hDPint->GetYaxis()->GetBinCenter(hDPint->GetNbinsY()-1), hDPint->GetZaxis()->GetBinCenter(hDPint->GetNbinsZ()-1));
    // Interpolate the distortion
    float phi=pos.Phi();
    if (phi<0) phi+=2*TMath::Pi();
    assert(checkBounds(phi, pos.Perp(), pos.Z(), hDPint));
    double dp = hDPint->Interpolate(phi, pos.Perp(), pos.Z());
    double dr = hDRint->Interpolate(phi, pos.Perp(), pos.Z());
    double dz = hDZint->Interpolate(phi, pos.Perp(), pos.Z());
    // To get the vector distortion, we need to find the delta between the original vector and the vector + distortion
    //printf("position=%f %f %f", pos.Phi(), pos.Perp(), pos.Z());
    //printf("  dp: %f dr: %f dz: %f\n", dp, dr, dz);
    if (!isRadians){
        dp/=pos.Perp();
    }
    TVector3 pfinal = pos;
    pfinal.SetPhi(pos.Phi() - dp);
    pfinal.SetPerp(pos.Perp() - dr);
    pfinal.SetZ(pos.Z() - dz);
    return pfinal;
}


void compositeCorrection(std::string firstfile, std::string secondfile){
    // Open the first correction in order:
    TFile *modulecorr_tfile = new TFile(firstfile.c_str());
    if (!modulecorr_tfile->IsOpen())
    {
        std::cout << "Error: could not open file " << firstfile << std::endl;
        return;
    }

    // Open the second correction in order:

    TFile *distortion_tfile = new TFile(secondfile.c_str());
    if (!distortion_tfile->IsOpen())
    {
        std::cout << "Error: could not open file " << secondfile << std::endl;
        return;
    }

    // get the module boundary correction:
        const std::array<const std::string, 2> extension = {{"_negz", "_posz"}};

    TH3* hDPmod[2], *hDRmod[2], *hDZmod[2];
 for (int j = 0; j < 2; ++j)
    {
        hDPmod[j] = dynamic_cast<TH3*>(distortion_tfile->Get((std::string("hIntDistortionP")+extension[j]).c_str()));
        assert(hDPmod[j]);
        hDRmod[j] = dynamic_cast<TH3*>(distortion_tfile->Get((std::string("hIntDistortionR")+extension[j]).c_str()));
        assert(hDRmod[j]);
        hDZmod[j] = dynamic_cast<TH3*>(distortion_tfile->Get((std::string("hIntDistortionZ")+extension[j]).c_str()));
        assert(hDZmod[j]);   
    }

    // Get the static correction:
    TH3* hDPint[2], *hDRint[2], *hDZint[2];
    for (int j = 0; j < 2; ++j)
    {
        hDPint[j] = dynamic_cast<TH3*>(distortion_tfile->Get((std::string("hIntDistortionP")+extension[j]).c_str()));
        assert(hDPint[j]);
        hDRint[j] = dynamic_cast<TH3*>(distortion_tfile->Get((std::string("hIntDistortionR")+extension[j]).c_str()));
        assert(hDRint[j]);
        hDZint[j] = dynamic_cast<TH3*>(distortion_tfile->Get((std::string("hIntDistortionZ")+extension[j]).c_str()));
        assert(hDZint[j]);   
    }


//copy the static correction histograms to new histograms in the new file, in preparation for the composite correction
    TFile *f = new TFile("compositeCorrection.static.and.module.root","recreate");
    TH3* hDPcomposite[2], *hDRcomposite[2], *hDZcomposite[2];

    for (int j = 0; j < 2; ++j)
    {
        hDPcomposite[j] = dynamic_cast<TH3*>(hDPint[j]->Clone());
        hDRcomposite[j] = dynamic_cast<TH3*>(hDRint[j]->Clone());
        hDZcomposite[j] = dynamic_cast<TH3*>(hDZint[j]->Clone());
    }

//then rename the original histograms for easy identification
    for (int j = 0; j < 2; ++j)
    {
        hDPint[j]->SetName((std::string("staticP")+extension[j]).c_str());
        hDRint[j]->SetName((std::string("staticR")+extension[j]).c_str());
        hDZint[j]->SetName((std::string("staticZ")+extension[j]).c_str());
        hDPmod[j]->SetName((std::string("moduleP")+extension[j]).c_str());
        hDRmod[j]->SetName((std::string("moduleR")+extension[j]).c_str());
        hDZmod[j]->SetName((std::string("moduleZ")+extension[j]).c_str());

    }


    //loop over the bin centers of the static correction histograms
    for (int i=1; i<=hDPcomposite[0]->GetNbinsY(); i++){
 
        float rpos=hDPcomposite[0]->GetYaxis()->GetBinCenter(i);
        if (i==1) rpos=hDPcomposite[0]->GetYaxis()->GetBinCenter(i+1);
        if (i==hDPcomposite[0]->GetNbinsY()) rpos=hDPcomposite[0]->GetYaxis()->GetBinCenter(i-1);
        TVector3 pos(rpos,0,0);
        for (int j=1; j<=hDPcomposite[0]->GetNbinsX(); j++){
            float phipos=hDPcomposite[0]->GetXaxis()->GetBinCenter(j);
            pos.SetPhi(phipos);
            for (int k=1; k<=hDPcomposite[0]->GetNbinsZ(); k++){
                float zpos=hDPcomposite[0]->GetZaxis()->GetBinCenter(k);
                if (k==1) zpos=hDRcomposite[0]->GetZaxis()->GetBinCenter(k+1);
                if (k==hDPcomposite[0]->GetNbinsZ()) zpos=hDPcomposite[0]->GetZaxis()->GetBinCenter(k-1);

                for (int side=0;side<2;side++){
                    if (side>0)zpos*=-1;//swap signs
                    pos.SetZ(zpos);
                    //get the module edge correction at the raw position
                    TVector3 pos1 = correctPosition(pos, hDPmod[side], hDRmod[side], hDZmod[side], true);
                    //get the static correction at the corrected position
                    TVector3 pos2 = correctPosition(pos1, hDPint[side], hDRint[side], hDZint[side], false);
                    //get the total correction, in terms of r, phi, and z
                    float dphi = pos2.Phi()-pos.Phi();//might have a twopi anomaly here...
                    float dr = pos2.Perp()-pos.Perp();
                    float dz = pos2.Z()-pos.Z();

                    //set the correction in the composite histograms
                    //nb this is in radians not cm units for phi.
                    hDPcomposite[side]->SetBinContent(j,i,k,-dphi);
                    hDRcomposite[side]->SetBinContent(j,i,k,-dr);
                    hDZcomposite[side]->SetBinContent(j,i,k,-dz);
                }
            }
        }
    }
//write the composite histograms to the new file
    for (int j = 0; j < 2; ++j)
    {
        hDPcomposite[j]->Write();
        hDRcomposite[j]->Write();
        hDZcomposite[j]->Write();
    }
    f->Close();


    return;
}