        "gamma_is_med",
        "gamma_is_fake",
        "pass_selection1",
        "pass_selection2",
        "photon_selection",
        "med_photon_pt",
        "med_photon_eta",
        "med_photon_phi",
        "med_photon_sieie",
        "med_photon_chiso",
        "med_photon_is_barrel",
        "med_photon_is_endcap",
        "med_photon_gen_matching",
        "med_photon_isprompt",
        "fake_photon_pt",
        "fake_photon_eta",
        "fake_photon_phi",
        "fake_photon_sieie",
        "fake_photon_chiso",
        "fake_photon_is_barrel",
        "fake_photon_is_endcap",
        "fake_photon_gen_matching",
        "fake_photon_isprompt",
    
    def select_fake_and_medium_photons(self, events, photons):
        # | pt | scEta | H over EM | sigma ieie | Isoch | IsoNeu | Isopho | 
        # listed from the right side
        mask1 = 0b10101010101010  # full medium ID
        mask2 = 0b00101010101010  # remove Isopho
        mask3 = 0b10001010101010  # remove IsoNeu
        mask4 = 0b10100010101010  # remove Isoch 
        mask5 = 0b10101000101010  # remove sigma ieie
        mask6 = 0b10100000101010  # remove the Isoch and sigma ieie

        # photons = photons[photons.pixelSeed]
        bitmap = photons.vidNestedWPBitmap

        # select medium and control photons
        # after adding the photons that pass the full ID, add the photons that pass the inverted ID
        # select control photons that don't pass the full ID but pass ID that one of cut inverted, which means this cut is inverted
        medium_and_control_cut = ((bitmap & mask1) == mask1) | ((bitmap & (mask2 + (3<<12))) == mask2) | ((bitmap & (mask3 + (3<<10))) == mask3) | ((bitmap & (mask4 + (3<<8))) == mask4) | ((bitmap & (mask5 + (3<<6))) == mask5)
        selected_medium_or_control_photons = photons[medium_and_control_cut] # append the medium and control photons

        # select fake photons
        # save photons pass the ID without sigma ieie and Isoch
        fake_template_photons_cut = (bitmap & (mask6 + (3<<6) + (3<<8))) == mask6
        selected_fake_template_photons = photons[fake_template_photons_cut]  #for fake template from data

        pass_selection1 = awkward.num(selected_medium_or_control_photons) >= 1  # select medium and control photons without fake photon
        pass_selection2 = awkward.num(selected_fake_template_photons) >= 1  # select fake photons

        # pass_selection1 and pass_selection2 can appear meantime
        awkward_utils.add_field(events, "pass_selection1", pass_selection1)  # has selected medium and control photons
        awkward_utils.add_field(events, "pass_selection2", pass_selection2)  # has selected fake photons

        isprompt_mask = (1 << 0)  #isPrompt
        isdirectprompttaudecayproduct_mask = (1 << 5)  #isDirectPromptTauDecayProduct
        isfromhardprocess_mask = (1 << 8)  #isPrompt
        isdirecttaudecayproduct_mask = (1 << 4)  #isDirectTauDecayProduct
        isprompttaudecayproduct = (1 << 3)  #isPromptTauDecayProduct

        med_cand = awkward.firsts(selected_medium_or_control_photons)

        med_photon_is_barrel = (med_cand.eta < self.options["photons"]["eta"][0][1]) & (med_cand.eta > self.options["photons"]["eta"][0][0])
        med_photon_is_endcap = (med_cand.eta < self.options["photons"]["eta"][1][1]) & (med_cand.eta > self.options["photons"]["eta"][1][0])
        awkward_utils.add_field(events, "med_photon_is_barrel", awkward.fill_none(med_photon_is_barrel, -1))
        awkward_utils.add_field(events, "med_photon_is_endcap", awkward.fill_none(med_photon_is_endcap, -1))

        for field in ["pt", "eta", "phi", "sieie"]:
            awkward_utils.add_field(
                events,
                "med_photon_%s" % field,
                awkward.fill_none(getattr(med_cand, field), DUMMY_VALUE)
            )
        awkward_utils.add_field(events, "med_photon_chiso",  awkward.fill_none(med_cand.pfRelIso03_chg, DUMMY_VALUE))

        if hasattr(events, 'nGenPart') and hasattr(med_cand,'genPartIdx'):
            genparts = events.GenPart
            match_gen = awkward.local_index(genparts, axis=1) == med_cand.genPartIdx
            mother_match_gen = awkward.local_index(genparts, axis=1) == awkward.sum(genparts[match_gen].genPartIdxMother, axis=1)

            gen_cut1 = (med_cand.genPartIdx >= 0) & (genparts[match_gen].pdgId  == 22)
            gen_cut2 = ((genparts[match_gen].statusFlags & isprompt_mask == isprompt_mask) | (genparts[match_gen].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)) & (genparts[match_gen].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)
            gen_cut3 = ((genparts[match_gen].statusFlags & isprompt_mask == isprompt_mask) | (genparts[match_gen].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask))
            gen_cut4 = (genparts[match_gen].genPartIdxMother >= 0 & (awkward.sum(abs(genparts[mother_match_gen].pdgId) == 11, axis=1) | awkward.sum(abs(genparts[mother_match_gen].pdgId) == 13, axis=1) | awkward.sum(abs(genparts[mother_match_gen].pdgId) == 15, axis=1)))
            gen_cut5 = (med_cand.genPartIdx >= 0) & (abs(genparts[match_gen].pdgId) == 11)

            photon_gen_matching = (gen_cut1 & gen_cut2) * 6 + (gen_cut1 & gen_cut3 & ~gen_cut4) * 5 + (gen_cut1 & gen_cut3 & gen_cut4) * 4 + (gen_cut1 & (~gen_cut2 & ~gen_cut3)) * 3 + (gen_cut5 & ~gen_cut3) * 2 + (gen_cut5 & gen_cut3) * 1 + (~gen_cut1 & ~gen_cut5) * 0

            photon_gen_matching = awkward.sum(photon_gen_matching, axis=1)
            awkward_utils.add_field(events, "med_photon_gen_matching",  awkward.fill_none(photon_gen_matching, -1))

            genpart_cut = (genparts.pt > 5) & (abs(genparts.pdgId) == 22) 
            photon_isprompt = (med_cand.genPartIdx >=0) &  ((genparts[match_gen].statusFlags & isprompt_mask == isprompt_mask) | (genparts[match_gen].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask) | (genparts[match_gen].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)) & (awkward.sum(~object_selections.delta_R(genparts[genpart_cut], med_cand, 0.3), axis=1) >= 1)

            photon_isprompt = awkward.sum(photon_isprompt, axis=1)
            awkward_utils.add_field(events, "med_photon_isprompt",  awkward.fill_none(photon_isprompt, -1))

        fake_cand = awkward.firsts(selected_fake_template_photons)

        fake_photon_is_barrel = (fake_cand.eta < self.options["photons"]["eta"][0][1]) & (fake_cand.eta > self.options["photons"]["eta"][0][0])
        fake_photon_is_endcap = (fake_cand.eta < self.options["photons"]["eta"][1][1]) & (fake_cand.eta > self.options["photons"]["eta"][1][0])
        awkward_utils.add_field(events, "fake_photon_is_barrel", awkward.fill_none(fake_photon_is_barrel, -1))
        awkward_utils.add_field(events, "fake_photon_is_endcap", awkward.fill_none(fake_photon_is_endcap, -1))

        for field in ["pt", "eta", "phi", "sieie"]:
            awkward_utils.add_field(
                events,
                "fake_photon_%s" % field,
                awkward.fill_none(getattr(fake_cand, field), DUMMY_VALUE)
            )
        awkward_utils.add_field(events, "fake_photon_chiso", awkward.fill_none(fake_cand.pfRelIso03_chg, DUMMY_VALUE))

        if hasattr(events, 'nGenPart') and hasattr(fake_cand,'genPartIdx'):
            genparts = events.GenPart
            match_gen = awkward.local_index(genparts, axis=1) == fake_cand.genPartIdx
            mother_match_gen = awkward.local_index(genparts, axis=1) == awkward.sum(genparts[match_gen].genPartIdxMother, axis=1)

            gen_cut1 = (fake_cand.genPartIdx >= 0) & (genparts[match_gen].pdgId  == 22)
            gen_cut2 = ((genparts[match_gen].statusFlags & isprompt_mask == isprompt_mask) | (genparts[match_gen].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask)) & (genparts[match_gen].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)
            gen_cut3 = ((genparts[match_gen].statusFlags & isprompt_mask == isprompt_mask) | (genparts[match_gen].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask))
            gen_cut4 = (genparts[match_gen].genPartIdxMother >= 0 & (awkward.sum(abs(genparts[mother_match_gen].pdgId) == 11, axis=1) | awkward.sum(abs(genparts[mother_match_gen].pdgId) == 13, axis=1) | awkward.sum(abs(genparts[mother_match_gen].pdgId) == 15, axis=1)))
            gen_cut5 = (fake_cand.genPartIdx >= 0) & (abs(genparts[match_gen].pdgId) == 11)

            photon_gen_matching = (gen_cut1 & gen_cut2) * 6 + (gen_cut1 & gen_cut3 & ~gen_cut4) * 5 + (gen_cut1 & gen_cut3 & gen_cut4) * 4 + (gen_cut1 & (~gen_cut2 & ~gen_cut3)) * 3 + (gen_cut5 & ~gen_cut3) * 2 + (gen_cut5 & gen_cut3) * 1 + (~gen_cut1 & ~gen_cut5) * 0

            photon_gen_matching = awkward.sum(photon_gen_matching, axis=1)
            awkward_utils.add_field(events, "fake_photon_gen_matching",  awkward.fill_none(photon_gen_matching, -1))

            genpart_cut = (genparts.pt > 5) & (abs(genparts.pdgId) == 22) 
            photon_isprompt = (fake_cand.genPartIdx >=0) &  ((genparts[match_gen].statusFlags & isprompt_mask == isprompt_mask) | (genparts[match_gen].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask) | (genparts[match_gen].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)) & (awkward.sum(~object_selections.delta_R(genparts[genpart_cut], fake_cand, 0.3), axis=1) >= 1)

            photon_isprompt = awkward.sum(photon_isprompt, axis=1)
            awkward_utils.add_field(events, "fake_photon_isprompt",  awkward.fill_none(photon_isprompt, -1))


        bitmap = med_cand.vidNestedWPBitmap 
        photon_selection = ((bitmap & mask1) == mask1) * 1 + ((bitmap & (mask2 + (3<<12))) == mask2) * 2 + ((bitmap & (mask3 + (3<<10))) == mask3) * 3 + ((bitmap & (mask4 + (3<<8))) == mask4) * 4 + ((bitmap & (mask5 + (3<<6))) == mask5) * 5
        awkward_utils.add_field(events, "photon_selection", awkward.fill_none(photon_selection, -1))
        #pass_selection1 && photon_selection==5 -> build ture template from MC and data template from data
        #pass_selection1 && ((photon_selection!=1 && photon_selection==2) || (!1 && ==3) || (!1 && ==4) || (!=1 && ==5)) -> build fake photon enriched sample from data



            mother_match_gen = awkward.local_index(genparts, axis=1) == awkward.sum(genparts[match_gen].genPartIdxMother, axis=1)

            gen_cut1 = (med_cand.genPartIdx >= 0) & (genparts[match_gen].pdgId  == 22)
            gen_cut2 = (genparts[match_gen].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)
            gen_cut3 = ((genparts[match_gen].statusFlags & isprompt_mask == isprompt_mask) | (genparts[match_gen].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask))
            gen_cut4 = (genparts[match_gen].genPartIdxMother >= 0 & (awkward.sum(abs(genparts[mother_match_gen].pdgId) == 11, axis=1) | awkward.sum(abs(genparts[mother_match_gen].pdgId) == 13, axis=1) | awkward.sum(abs(genparts[mother_match_gen].pdgId) == 15, axis=1)))
            gen_cut5 = (med_cand.genPartIdx >= 0) & (abs(genparts[match_gen].pdgId) == 11)

            photon_gen_matching = (gen_cut1 & ~gen_cut2 & gen_cut3) * 7 + (gen_cut1 & gen_cut2 & gen_cut3) * 6 + (gen_cut1 & gen_cut3 & ~gen_cut4) * 5 + (gen_cut1 & gen_cut3 & gen_cut4) * 4 + (gen_cut1 & ~gen_cut3) * 3 + (gen_cut5 & ~gen_cut3) * 2 + (gen_cut5 & gen_cut3) * 1

            photon_gen_matching = awkward.sum(photon_gen_matching, axis=1)
            awkward_utils.add_field(events, "photon_gen_matching",  awkward.fill_none(photon_gen_matching, -1))


# | pt | scEta | H over EM | sigma ieie | Isoch | IsoNeu | Isopho | 
# listed from the right side
mask1 = 0b10101010101010  # full medium ID
mask2 = 0b00101010101010  # remove Isopho
mask3 = 0b10001010101010  # remove IsoNeu
mask4 = 0b10100010101010  # remove Isoch 
mask5 = 0b10101000101010  # remove sigma ieie
mask6 = 0b10100000101010  # remove the Isoch and sigma ieie

bitmap = photons.vidNestedWPBitmap

medium_and_control_cut = ((bitmap & mask1) == mask1) | ((bitmap & (mask2 + (3<<12))) == mask2) | ((bitmap & (mask3 + (3<<10))) == mask3) | ((bitmap & (mask4 + (3<<8))) == mask4) | ((bitmap & mask6) == mask6)
selected_medium_or_control_photons = photons[medium_and_control_cut]
med_cand = awkward.firsts(selected_medium_or_control_photons)

bitmap = med_cand.vidNestedWPBitmap 
photon_selection = (
    ((bitmap & mask1) == mask1) * 1 + 
    ((bitmap & (mask2 + (3<<12))) == mask2) * 2 + 
    ((bitmap & (mask3 + (3<<10))) == mask3) * 3 + 
    ((bitmap & (mask4 + (3<<8))) == mask4) * 4 + 
    ((bitmap & (mask5 + (3<<6))) == mask5) * 5 + 
    ((bitmap & (mask6 + (3<<6) + (3<<8))) == mask6) * 6 + 
    ((((bitmap & mask5)) == mask5) & ((bitmap & mask1) != mask1) & ((bitmap & (mask5 + (3<<6))) != mask5)) * 7 + 
    ((bitmap & mask6 == mask6) & (((bitmap & mask5)) != mask5) & ((bitmap & (mask4 + (3<<8))) != mask4) & ((bitmap & (mask6 + (3<<6) + (3<<8))) != mask6)) * 8
)
awkward_utils.add_field(events, "photon_selection", awkward.fill_none(photon_selection, -1))

isprompt_mask = (1 << 0)  #isPrompt
isdirectprompttaudecayproduct_mask = (1 << 5)  #isDirectPromptTauDecayProduct
isfromhardprocess_mask = (1 << 8)  #isPrompt

if hasattr(events, 'nGenPart') and hasattr(med_cand,'genPartIdx'):
    match_gen = awkward.local_index(events.GenPart, axis=1) == med_cand.genPartIdx

    genpart_cut = (events.GenPart.pt > 5) & (abs(events.GenPart.pdgId) == 22) 
    photon_isprompt = (med_cand.genPartIdx >=0) &  ((events.GenPart[match_gen].statusFlags & isprompt_mask == isprompt_mask) | (events.GenPart[match_gen].statusFlags & isdirectprompttaudecayproduct_mask == isdirectprompttaudecayproduct_mask) | (events.GenPart[match_gen].statusFlags & isfromhardprocess_mask == isfromhardprocess_mask)) & (awkward.sum(~object_selections.delta_R(events.GenPart[genpart_cut], med_cand, 0.3), axis=1) >= 1)

    photon_isprompt = awkward.sum(photon_isprompt, axis=1)
    awkward_utils.add_field(events, "photon_isprompt",  awkward.fill_none(photon_isprompt, -1))