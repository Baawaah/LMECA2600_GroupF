# Auxiliary Function
# Code by the group F
# by Pierre-Yves Legros
#    Lise Ceresiat
#    Son Tran




function M = molarMass(X)
    X = [U235 U238 U239 Np239 Pu239 Xe135]

function demi_vie = Demi_vie (X,Transfo)
    X = [U235 U238 U239 Np239 Pu239 Xe135]
    Transfo = [Alpha,BetaPlus,BetaMinus]

function [sigma] = Section_efficace(X,Transfo,n_eV,User_Adress)
    X = [U235 U238 U239 Np239 Pu239]
    Transfo = [Fission Capture]
    n_eV = [eV]
    User_Adress [location of the files, e.g : data]
