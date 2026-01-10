let randomNumbers;

document.addEventListener("DOMContentLoaded", function () {
  const button = document.getElementById("generate-btn");
  const display = document.getElementById("random-number");

  // Generate the initial random number
  randomNumbers = Math.floor(Math.random() * 20 + 1);

  button.addEventListener("click", function () {
    randomNumbers = Math.floor(Math.random() * 20 + 1); // this will Update random number on click
    console.log(randomNumbers);
    //display.textContent = `Generated number is: ${randomNumbers}`; we ont wwant the teacher to know the answer..
  });
});

function checkValue() {
  let userGuess = $(".form-control").val();
  userGuess = parseInt(userGuess);
   let guesses = [];
  guesses.push(userGuess);
  $(".guessArray").append(userGuess + ", ").css("border", "orange solid 3px");
  

  console.log("before parseInt(): " + typeof userGuess);
  userGuess = parseInt(userGuess);
  console.log("after parseInt(): " + typeof userGuess);

  if (userGuess === randomNumbers) {
    $("#result").text(`You guessed it right! It took you${guesses.length} tries`).css("color", "green");
  } else if (userGuess < randomNumbers) {
    $("#result").text("Number you picked is too low").css("color", "red");
  } else if (userGuess > randomNumbers) {
    $("#result").text("Number you picked is too high!").css("color", "red");
  } else {
    $("#result").text("Please enter a number!").css("color", "red");
  }
}
