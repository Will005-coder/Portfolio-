#include <ESP32Encoder.h>
#include <esp_timer.h>
// #include //wtv the jetson is telling

#define IN1 6
#define IN2 7

#define PWM 5

#define CHANNEL_A 11
#define CHANNEL_B 12

ESP32Encoder encoder;

//Encoder specification
const float ENCODER_PPR = 751.8f;
const int QUADRATURE_MULTIPLIER = 4; //We're cointing 4 edges/signals per click on the encoder(higher-res)
const long COUNTS_PER_REVOLUTION = (long)(ENCODER_PPR * QUADRATURE_MULTIPLIER);

//How often we read the encoder (100 ms)
const uint64_t READ_INTERVAL_MICROSECONDS = 100000ULL;

//These values change spontaneously
volatile int64_t latestEncoderCount = 0;
volatile int64_t latestTimestampMicroseconds = 0;
volatile bool newReadingAvailable = false;

//Previous reading (used to compute change)
int64_t prevEncoderCount = 0;
int64_t prevTimestampMicroseconds = 0;

float prevRPM = 0;
float currentRPMWeight = 0;
float cleanRPM = 0;

double dt, last_time;
double integral, previousError, pidOutput = 0; 
//PID params
double Kp, Ki, Kd; //for now just following a setup tutorial, we can tune this later
double setpoint = 10;

//Periodic-timer Callback
static void IRAM_ATTR timerCallback(void* arg) {
  latestTimestampMicroseconds = esp_timer_get_time();
  latestEncoderCount = encoder.getCount();
  newReadingAvailable = true;
}


// ---------- make RPM calc a function ----------
double calculateRpmAndFilter() {
  //Compute changes in encoder count and time (using globals set by timer)
  int64_t countChange = latestEncoderCount - prevEncoderCount;
  int64_t timeChangeMicroseconds = latestTimestampMicroseconds - prevTimestampMicroseconds;

  // protect against zero/negative dt
  if (timeChangeMicroseconds <= 0) {
    return prevRPM; // nothing new, return previous filtered RPM
  }

  // Correct RPM formula using microseconds -> minutes:
  // RPM = (countChange * 60,000,000) / (countsPerRev * deltaTimeUs)
  double numerator = (double)countChange * 60000000.0;
  double denominator = (double)COUNTS_PER_REVOLUTION * (double)timeChangeMicroseconds;
  double rawRPM = 0.0;
  if (denominator != 0.0) rawRPM = numerator / denominator;

  // Exponential Moving Average: currentRPM_Weight used to smooth
  // currentRPMWeight should be 0..1 (0 = max smoothing, 1 = no smoothing)
  cleanRPM = (currentRPMWeight * (float)rawRPM) + ((1.0f - currentRPMWeight) * prevRPM);

  // store raw->prev
  prevRPM = cleanRPM;

  return cleanRPM;
}


void setup() {
  //put your setup code here, to run once:
  Kp = 0;
  Ki = 0;
  Kd = 0;
  Serial.begin(9600);

  pinMode(PWM,OUTPUT);
  pinMode(IN1,OUTPUT);
  pinMode(IN2,OUTPUT);
  //idk random rpm right now
  setpoint = 10;
  for(int i = 0; i < 50; i++)
  {
    Serial.print(setpoint);
    Serial.print(",");
    Serial.println(0);
  }
  Serial.begin(115200);
  Serial.println("Encoder RPM - cleaned variable names");

  pinMode(CHANNEL_A, INPUT_PULLUP);
  pinMode(CHANNEL_B, INPUT_PULLUP);

  encoder.attachFullQuad(CHANNEL_A, CHANNEL_B);
  encoder.clearCount();
  delay(10);

  prevTimestampMicroseconds = esp_timer_get_time();
  prevEncoderCount = encoder.getCount();

  //Create periodic timer that will run on the esp32, using the callback from timerCallback
  const esp_timer_create_args_t timerArgs = {
    .callback = &timerCallback,
    .name = "encoder_reader"
  };

  esp_timer_handle_t periodicTimer;
  esp_timer_create(&timerArgs, &periodicTimer);
  esp_timer_start_periodic(periodicTimer, READ_INTERVAL_MICROSECONDS);

  // initialize PID timing
  last_time = (double) millis();
  integral = 0.0;
  previousError = 0.0;

  // set a reasonable EMA weight default (you can tune)
  currentRPMWeight = 0.2f;
}



void loop() {
  if (Serial.available() > 0) {
    double incoming = Serial.parseFloat();  // read number user typed
    if (incoming > 0) {                     // simple safety: ignore zeros
      setpoint = incoming;
      Serial.print("New setpoint -> ");
      Serial.println(setpoint);
    }
  }

  if (!newReadingAvailable) return;
  
  //Pause the rest of the code to refresh start and end encodercount and seconds
  noInterrupts();
  int64_t currentEncoderCount = latestEncoderCount;
  int64_t currentTimestampMicroseconds = latestTimestampMicroseconds;
  newReadingAvailable = false;
  interrupts();

  //Compute changes in encoder count and time
  int64_t countChange = currentEncoderCount - prevEncoderCount;
  int64_t timeChangeMicroseconds = currentTimestampMicroseconds - prevTimestampMicroseconds;

  //These tell us if the encoder is even counting
  Serial.print("deltacount: ");
  Serial.print(countChange);
  Serial.print(" | deltaT(us): ");
  Serial.println(timeChangeMicroseconds);

  //Change previousEncoderCount refernce point once time passes
  if (timeChangeMicroseconds <= 0) {
    prevEncoderCount = currentEncoderCount;
    prevTimestampMicroseconds = currentTimestampMicroseconds;
    return;
  }

  //Integer RPM calculation:
  //RPM = (countChange * 6,000,000) / (COUNTS_PER_REVOLUTION * timeChangeMicroseconds)
  int64_t numerator = countChange * 60000000LL;  
  int64_t denominator = (int64_t)COUNTS_PER_REVOLUTION * timeChangeMicroseconds;

  bool isNegative = false;
  if (numerator < 0) {
    isNegative = true;
    numerator = -numerator;
  }

  //Compute RPM*10 (one decimal place)
  int64_t rpmTimes10 = (numerator * 10 + denominator / 2) / denominator;

  if (isNegative) rpmTimes10 = -rpmTimes10;

  //Exponential Moving Average: currentRPM_Weight = 0.4
  currentRPMWeight = 0.2;

  //clean our noise with EMA
  float rawRPM = rpmTimes10 / 10.0f;
  cleanRPM = (rawRPM * currentRPMWeight) + (prevRPM * (1 - currentRPMWeight));

  // DEBUG PRINT
  Serial.print("Raw RPM: ");
  Serial.print(rawRPM);
  Serial.print(" | Clean RPM: ");
  Serial.println(cleanRPM);

  Serial.print("RPM: ");
  Serial.print(rpmTimes10 / 10);
  Serial.print('.');
  Serial.println(llabs(rpmTimes10 % 10));

  //Store for next reading
  prevEncoderCount = currentEncoderCount;
  prevTimestampMicroseconds = currentTimestampMicroseconds;
  prevRPM = rawRPM;

 //Timer keeping track of time elapsed.  
  double now = millis();
  dt = (now-last_time)/1000.00;
  last_time = now;

  // float setpoint = //we need to request a setpoint

  //Calculates the error in the system
  double actualRpm = calculateRpmAndFilter(); // from our rpm calc
  double error = setpoint - actualRpm;
  pidOutput = pid(error, dt);

  //Where voltage and direction get set
  int dir = 1;
  if (pidOutput < 0){ dir=-1; }
  //pwr is the new voltage, which is the absolute value of the output. 
  int pwr = (int) fabs(pidOutput);
  if(pwr > 255) pwr = 255;
  setMotor(dir,pwr,PWM,IN1,IN2);

  //print out values
  Serial.print("Setpoint: "); Serial.print(setpoint);
  Serial.print(", Actual: "); Serial.print(actualRpm);
  Serial.print(", Error: "); Serial.println(error);

}

// note: kept your pid function name and comments, fixed internals
double pid(double error, double dt){
  //PID calculation
  //First gets the proportional part of PID, which is just the error.
  double proportional = error;

  // The integral part sums up all the error and multiplies it by change in time. (Rough approx)
  // add clamping
  integral += error * dt;
  const double INTEGRAL_LIMIT = 1000.0;
  //Set integral limit: CAN'T GO PAST 1000RPM
  if (integral > INTEGRAL_LIMIT) integral = INTEGRAL_LIMIT;
  if (integral < -INTEGRAL_LIMIT) integral = -INTEGRAL_LIMIT;

  //derivative finds the slope of the errors (Rough approx)
  double rawDerivative = 0.0;
  if (dt > 0.0) rawDerivative = (error - previousError) / dt;

  //creates the variable derivativeLPF which is used as the refrence point in LPF calculation
  //set equal to 0 only once. 
  static double derivativeLPF = 0;
  //weight is what effects how much the real derivative effects the value we pass to the PID. 0.1 is used because it offers a smooth but not too smooth value. 
  double weight = 0.1; //this val can be tuned
  //LPF is the preivous lpf times the difference in previous derivative and current derivative multiplied by a small amount. 
  derivativeLPF = derivativeLPF + weight * (rawDerivative - derivativeLPF);
  double derivative  = derivativeLPF;
  previousError = error;
  //Calculates new output
  double output = (Kp * proportional) + (Ki * integral) + (Kd * derivative);
  return output;
}

//function that takes the direction, the pwrm, and the two motor pins and changes the speed of motor accordingly. 
void setMotor(int dir, int pwmVal, int pwm, int in1, int in2){
  analogWrite(pwm,pwmVal);
  if(dir==1){ digitalWrite(in1,HIGH); digitalWrite(in2,LOW); }
  else if(dir==-1){ digitalWrite(in1,LOW); digitalWrite(in2,HIGH); }
  else{ digitalWrite(in1,LOW); digitalWrite(in2,LOW); }
}
