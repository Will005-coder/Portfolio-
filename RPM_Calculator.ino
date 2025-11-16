#include <ESP32Encoder.h>
#include <esp_timer.h>

#define CHANNEL_A 11
#define CHANNEL_B 12

ESP32Encoder encoder;

// Set this to your encoder spec (pulses per channel per revolution)
const float PPR = 751.8;          // Number of Clicks/pulses per rotation from manufaturers sheet
const int QUAD_MULT = 4;          // 2 edges per pulse for attachHalfQuad, 4 edges per pulse for attachFullQuad (We're using this for higher resolution)
const long COUNTS_PER_REV = (PPR * QUAD_MULT);

int64_t lastCount = 0;            // Tracks encoder count from last measurement
int64_t lastTime = 0;             // Time since RPM was last recorded (use 64-bit for microseconds)

void setup() {
  Serial.begin(115200);
  Serial.println("Encoder rotation + RPM demo");

  pinMode(CHANNEL_A, INPUT_PULLUP);
  pinMode(CHANNEL_B, INPUT_PULLUP);

  encoder.attachFullQuad(CHANNEL_A, CHANNEL_B);
  // If you switch to half quad, change QUAD_MULT to 2 and call:
  // encoder.attachHalfQuad(CHANNEL_A, CHANNEL_B);

  encoder.clearCount();           // start clean
  lastTime = esp_timer_get_time();
  lastCount = encoder.getCount();
}

void loop() {
  // Read current encoder count
  int64_t currentMicros = esp_timer_get_time();
  int64_t count = encoder.getCount();

  // Calculate change in counts and time
  int64_t deltaCountsPer100Ms = count - lastCount;                 // Change in encoder counts
  int64_t deltaTime = currentMicros - lastTime;            // Change in time for RPM calculation (microseconds)

  float RPM = 0.0; // Reset RPM every loop for recalculation

  // RPM = (deltaCountsPer100Ms * 60 * 1e6) / (COUNTS_PER_REV * deltaTime)
  RPM = (float)deltaCountsPer100Ms / (float)COUNTS_PER_REV;

  // Only print every 100 ms (100,000 microseconds)
  if (deltaTime >= 6000000) { // We'll need to use a callback function if we decide to use the periodic timer
    Serial.print("RAW COUNT: ");
    Serial.println(count);
    Serial.print("RPM: ");
    Serial.println(RPM, 3);  // show 3 decimal places
    Serial.print("Delta Counts: ");
    Serial.println(deltaCountsPer100Ms);
    Serial.print("Delta Time (us): ");
    Serial.println((long)deltaTime);

    encoder.clearCount();           // start clean
    lastCount = count;         // Reset count reference
  }
}
