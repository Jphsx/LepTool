
main: E_bank_manager.C E_bank.h E_bank_fit.h
	g++ -o banktest E_bank_manager.C `root-config --cflags --libs`

