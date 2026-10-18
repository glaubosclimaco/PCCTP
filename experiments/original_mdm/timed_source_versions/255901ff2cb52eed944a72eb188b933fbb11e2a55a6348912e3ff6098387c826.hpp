#pragma once
void research_initialize(const char*,int,int,int,int,const char*);
void research_snapshot(int,double,int*);
void research_mining();
void research_export(const char*,double,int*);
void research_move(int,double,double,double,int*,int*);
void research_process_start();
bool research_budget_active();
bool research_expired();
void research_progress(double,int*);
