/*
READ ME
Script Name : ttMotorController.ino
This script is designed for the Acebott ESP32 MAX dev board.
It uses the Bluepad library to remotely control four TT motors through an 
L298N motor driver, handling motor direction control.

DEPENDENCIES:
- esp32 library V 3.0.0 by espressif  
- Bluepad32 library
- Copy and paste for in File > Preferences > Additional board manager URLS {
      https://arduino.esp8266.com/stable/package_esp8266com_index.json
      https://espressif.github.io/arduino-esp32/package_esp32_index.json
      https://raw.githubusercontent.com/espressif/arduino-esp32/gh-pages/package_esp32_index.json
      https://raw.githubusercontent.com/ricardoquesada/esp32-arduino-lib-builder/master/bluepad32_files/package_esp32_bluepad32_index.json
      https://www.arduino.me/package_esp32_index.json
      }

AUTHOR NOTES:
- Use find "@" (`Ctrl` + F ) to navigate between major code segments. 
- Set Serial monitor baud rate to 921600, or adjust in @Setup as necessary

Segment Headers are listed below
- @Motor pins | Acebott Pin assignement and variable naming conventions
- @main |  Main script function call (acts as script setup)
- @Setup | Acebott Pin and controller setups
- @Loop | Guard clauses for identifying live gamepads
- @Callbacks | Updates a list of active controllers upon connection and disconnection
- @Helper function | contains `stopAllMotors` which is used as an E-stop mechanism
- @Gamepad input processing | Translates gamepad input into directional control

IMPROVEMENTS:
 - Fail-safe for uncontrollable motor state:
      A rotary encoder and PID should be implemented to prevent motor overload in the 
      prescence of inhibitory forces to the motor-shaft rotation.
*/
//@Motor pins
//F indicates front wheels
int FIN1 = 23; //Front-right wheel IN
int FIN2 = 19; //Front-right wheel OUT
int FIN3 = 14; //Front-left wheel IN
int FIN4 = 25; //Front-left wheel OUT

//B indicates back wheels
int BIN1 = 18; //Back-left wheel IN
int BIN2 = 17; //Back-left wheel OUT
int BIN3 = 26; //Back-left wheel IN
int BIN4 = 27; //Back-right wheel OUT

//@main
//Controller pointer
ControllerPtr myControllers[BP32_MAX_CONTROLLERS];

//Forward declare
void processGamepad(ControllerPtr gamepad);
void stopAllMotors();

//@Setup
void setup() {
  Serial.begin(921600);

  //Motor pins
  pinMode(FIN1, OUTPUT); pinMode(FIN2, OUTPUT);
  pinMode(FIN3, OUTPUT); pinMode(FIN4, OUTPUT);
  pinMode(BIN1, OUTPUT); pinMode(BIN2, OUTPUT);
  pinMode(BIN3, OUTPUT); pinMode(BIN4, OUTPUT);

  //Bluepad32 setup
  BP32.setup(&onConnectedGamepad, &onDisconnectedGamepad);
  BP32.forgetBluetoothKeys();
}

//@Loop
void loop() {
  BP32.update();

  for (int i=0; i<BP32_MAX_CONTROLLERS; i++) { //BP32_MAX_CONTROLLERS is typically 4
    ControllerPtr ctl = myControllers[i]; //Appends possible controllers to a list

    
    Checks if a controller is:
    - In `myControllers` or callable;
    - A gamepad; 
    - And still connected to the Acebott ESP32

    Sets `ctl` (controller) as the target for input in function call 
    
    if (ctl && ctl->isConnected() && ctl->isGamepad()) {
      processGamepad(ctl);
    }
  }

  delay(10); //optimal delay time for quick input acquisition with reduced overload on serial
}

//@Callbacks
void onConnectedGamepad(ControllerPtr ctl) {
  for (int i=0; i<BP32_MAX_CONTROLLERS; i++) {
    if (myControllers[i] == nullptr) {//Appends new controllers upon connection to  `myControllers` 
      myControllers[i] = ctl;
      Serial.print("Controller connected at index ");
      Serial.println(i);
      break;
    }
  }
}

void onDisconnectedGamepad(ControllerPtr ctl) {
  for (int i=0; i<BP32_MAX_CONTROLLERS; i++) {
    if (myControllers[i] == ctl) { //Updates `myControllers` by searching for a controller that disconnects and deletes it 
      myControllers[i] = nullptr;
      Serial.print("Controller disconnected from index ");
      Serial.println(i);
      break;
    }
  }
}

//@Helper function
void stopAllMotors() {
  digitalWrite(FIN1, LOW); digitalWrite(FIN2, LOW);
  digitalWrite(FIN3, LOW); digitalWrite(FIN4, LOW);
  digitalWrite(BIN1, LOW); digitalWrite(BIN2, LOW);
  digitalWrite(BIN3, LOW); digitalWrite(BIN4, LOW);
}

//@Gamepad input processing
void processGamepad(ControllerPtr gamepad) {
  long raw_ly = gamepad->axisY();      // LY for forward/reverse
  long raw_rx = gamepad->axisRX();     // RX for left/right
                                       
  float fy = -(float)raw_ly / 512.0f;  // forward `fy` units | negative is backwards
  float fx = (float)raw_rx / 512.0f;   // right `fx` units | negative is left

  const float deadzone = 0.18f; //Equivalent of rest point on controller joystick

  //LEFT STICK = FORWARD / REVERSE
  if (fy > deadzone) {
    // Forward
    digitalWrite(FIN1,HIGH); digitalWrite(FIN2,LOW);
    digitalWrite(FIN3,HIGH); digitalWrite(FIN4,LOW);
    digitalWrite(BIN1,LOW);  digitalWrite(BIN2,HIGH);  //B-side reversed
    digitalWrite(BIN3,LOW);  digitalWrite(BIN4,HIGH);
    return;   //override turning 
  }
  else if (fy < -deadzone) {
    //Reverse
    digitalWrite(FIN1,LOW);  digitalWrite(FIN2,HIGH);
    digitalWrite(FIN3,LOW);  digitalWrite(FIN4,HIGH);
    digitalWrite(BIN1,HIGH); digitalWrite(BIN2,LOW);   //B-side reversed
    digitalWrite(BIN3,HIGH); digitalWrite(BIN4,LOW);
    return;
  }

  //RIGHT STICK = LEFT / RIGHT TURN
  if (fx > deadzone) {
    //right (FIN3 & 4-> Backwards, BIN1 & 2-> Forward)
    digitalWrite(FIN1,LOW);  digitalWrite(FIN2,HIGH);
    digitalWrite(FIN3,HIGH); digitalWrite(FIN4,LOW);
    digitalWrite(BIN1,LOW);  digitalWrite(BIN2,HIGH);   //B-side reversed
    digitalWrite(BIN3,HIGH);  digitalWrite(BIN4,LOW);
  }
  else if (fx < -deadzone) {
    //left
    digitalWrite(FIN1,HIGH);  digitalWrite(FIN2,LOW);
    digitalWrite(FIN3,LOW); digitalWrite(FIN4,HIGH);
    digitalWrite(BIN1,HIGH);  digitalWrite(BIN2,LOW);   //B-side reversed
    digitalWrite(BIN3,LOW);  digitalWrite(BIN4,HIGH);
  }
  else {
    //Neutral\default mode: stop
    digitalWrite(FIN1,LOW); digitalWrite(FIN2,LOW);
    digitalWrite(FIN3,LOW); digitalWrite(FIN4,LOW);
    digitalWrite(BIN1,LOW); digitalWrite(BIN2,LOW);
    digitalWrite(BIN3,LOW); digitalWrite(BIN4,LOW);
  }

  //Stop all motors manually with `A` button
  if (gamepad->a()) stopAllMotors();
}
