################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Utils.cpp \
../genius.cpp \
../grasp.cpp \
../instances.cpp \
../instancia_loader.cpp \
../mineracao.cpp \
../principal.cpp \
../random_provider.cpp \
../regras_reducao.cpp \
../teste.cpp 

C_SRCS += \
../mt19937ar.c 

OBJS += \
./Utils.o \
./genius.o \
./grasp.o \
./instances.o \
./instancia_loader.o \
./mineracao.o \
./mt19937ar.o \
./principal.o \
./random_provider.o \
./regras_reducao.o \
./teste.o 

CPP_DEPS += \
./Utils.d \
./genius.d \
./grasp.d \
./instances.d \
./instancia_loader.d \
./mineracao.d \
./principal.d \
./random_provider.d \
./regras_reducao.d \
./teste.d 

C_DEPS += \
./mt19937ar.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


