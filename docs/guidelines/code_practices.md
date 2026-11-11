No argparse unless necessary when code needs to be called as subprocess (relaxed constraint): configuration belongs in a typed config object, not in CLI flags. Variables are defined at the top of the document when part of the invariants, or after the if __name__ == "__main__" clause. If the number of parameters exceeds 10, we create a module config.py containing a dataclass Config, which is imported and passed in main.

The reason is reproducibility, not taste: a config object is serialized next to the results it produced (config.json, job.json), so a result directory says exactly what made it. A flag like --phi-lr 1e-5 is a second source of truth that disappears with the terminal scrollback.

Two carve-outs:
- A process boundary may take a single positional path to a config file — the worker entrypoint is invoked as `python -m deepkernelgp.worker.entrypoint <job.json>`. One positional argument does not need argparse; if such a boundary ever needs three, revisit this rule rather than hand-rolling a parser.
- Third-party scripts we call but do not own (MoReHub's setup.py and embeddings.py use argparse) are out of scope. We pass them arguments; we do not restyle them.


## Following the recipe



Code should be thoroughly documented, with clear comments and docstrings.


When possible, definitions of variables and data can be preceded by a format similar to this one:
```
# DD. POS_X
# pos_x = int
# interp. the X coordinate position of an object in a 2D space relative to the center of the object
pos_x0 = 0 #stub
pos_x1 = 1

# DD. POS_Y
# pos_y = int
# interp. the Y coordinate position of an object in a 2D space relative to the center of the object
pos_y0 = 0 #stub
pos_y1 = 1

# DD. VEL_X
# vel_x = float
# interp. the x velocity of an object in a 2D space as a normalized vector
# Invariant: float within range -1,1
vel_x0 = 0 #stub
vel_x1 = 0.5

# DD. VEL_Y
# vel_y = float
# interp. the y velocity of an object in a 2D space as a normalized vector
# Invariant: float within range -1,1
vel_y0 = 0 #stub
vel_y1 = 0.5

# CDD. BALL
# ball = Ball(POS_X, POS_Y, VEL_X, VEL_Y)
# interp. a ball rendered in a 2D space with a position in x,y and velocity dx and dy
class Ball:
    def __init__(self,pos_x, pos_y, vel_x, vel_y):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.vel_x = vel_x
        self.vel_y = vel_y
    
    # FD simple_bounce
    # Signature: BALL.x, BALL.y -> None
    def simple_bounce(self):
        '''Test simulation of a ball bouncing in the x and y axis
        Algorithm:
            Take BALL.x and BALL.y and multiply them by -1
        '''
        self.pos_x *= self.vel_x
        self.pos_y *= self.vel_y

ball0 = Ball(pos_x0, pos_y0, vel_x0, vel_y0) #stub: base case
ball1 = Ball(pos_x1, pos_y1, vel_x1, vel_y1)

# DD. BALLS
# balls = [BALL, ..., n]
# interp. a collection of balls, where index reflect the order of appearance
balls0 = []
balls1 = [ball0]
balls2 = [ball0,ball1]


# DD. SCORE
# score = int
# interp. a natural number tracks the number of blocks destroyed in the game
# Invariant: SCORE is Natural, bigger or equal than 0
score0 = 0
score1 = 1

# DD. IS_GAME_RUNNING
# is_game_running = bool
# interp. Whether the game is running or not
is_game_running_0 = False
is_game_running_1 = True



# CDD. GAME
# game = Game()
# interp. an encapsulation of an entire game of 1976's Breakout game

class Game:
    def __init__(self):
        self.balls = [ball0, ball1]
        self.score = 0
        self.is_game_running = False
    
    # FD start_game
    # Signature: None -> None
    def start_game(self):
        '''Begin the game by placing the balls and setting IS_GAME_RUNNING to True
        Algorithm:
            for each BALL in BALLS:
                draw the BALL
            set IS_GAME_RUNNING to True
        '''
        for ball in self.balls:
            self.draw(ball)
        self.is_game_running = True

    # FD draw
    # Signature: BALL -> None
    def draw(self, obj):
        '''Draw the components making the game of Breakout
        Algorithm:
            draw each ball in balls
        '''
        obj.draw()

    #FD. update_balls():
    # Signature: None -> None
    def update_balls(self):
        '''update the position of the balls'''
        for ball in self.balls:
            ball.simple_bounce()


game0 = Game()
game1 = Game()

# DD. GAMES
# games = [GAME, ..., m]
# interp. a collection of games running simultaneously
games = []
games1 = [game0]
games2 = [game0, game1]

```

Programs should be built from clear templates that naturally compound into complex systems, allowing for clear separation of concerns and scalability. Each module produced this way should then be clearly tested. The implementation begins with a clear definition on the information domain: what motifs of information enter the program? We begin with atomic units of information, which further aggregate into arrays, dictionaries, class instances, etc. Each of these objects has a template to work with. Let's look at a couple of examples using the definitions previously outlined:


### Template for class 

```
# Template for BALL:
... (ball)
... ball.pos_x
... ball.pos_y
... ball.vel_x
... ball.vel_y

```

Here there are two main way to access the components of the ball. We can pass the ball as an object without accessing its attributes, where ... is replaced by a defined function, followed by the arguments of that function, in which we find ball:

```
... (ball)
```

This is the desirable behavior, which should promote the use of getters and setters to ensure proper encapsulation and avoid information leakage.

Interestingly, this template also gives potential universal templates in the call of functions

```

def fn_for_ball_like_objects(obj, *args):
    ... obj.draw()

fn_for_ball_like_objects(ball,*args)
fn_for_ball_like_objects(paddle,*args)

```
Moreover, it can give use a tool to execute functions within a universal template

### Template for higher-order functions

A higher-order function takes another function as an argument. It earns its place
only when it does something the passed function cannot do on its own, typically:

1. **It owns the traversal** of a collection, so the loop is written once.
2. **It enforces shared rules** (invariants, logging, bookkeeping) after every
   operation, so no individual function can forget them.

A wrapper that only forwards its arguments, such as `fn(ball, *args)`, adds a name
without adding behaviour. Call the function directly instead.

```
# Template for higher-order function over BALLS
def fn_for_balls(balls, fn, *args):
    for ball in balls:
        ... fn(ball, *args)
        ... (ball)
```

Any function passed as `fn` must follow the contract `BALL, ... -> None`. It takes
a BALL as its first argument, then its own arguments, and modifies the BALL in place.

Example:

```
# FD set_pos
# Signature: BALL, POS_X, POS_Y -> None
def set_pos(ball, pos_x, pos_y):
    '''Place BALL at (POS_X, POS_Y).
    Algorithm:
        delegate to BALL's position setter with the added functionality of checking invariants
    '''
    if pos_x >=0 and pos_y >=0:
        ball.set_pos(pos_x, pos_y)
    else:
        raise ValueError("The position of the ball has to be positive!)

# FD set_vel
# Signature: BALL, VEL_X, VEL_Y -> None
def set_vel(ball, vel_x, vel_y):
    '''Give BALL the velocity (VEL_X, VEL_Y).
    Algorithm:
        delegate to BALL's velocity setter
    '''
    if ball.pos_x + vel_x > SCREEN["X"] or ball.pos_x + vel_x < 0:
        ball.set_vel(-vel_x, vel_y)
    if ball.pos_y + vel_y > SCREEN["Y"]:
        ball.set_vel(vel_x, -vel_y)

# FD fn_for_balls
# Signature: BALLS, (BALL, ... -> None), ... -> None
def fn_for_balls(balls, fn, *args):
    '''Apply fn to every BALL in BALLS, then verify that each BALL
    still satisfies its invariants.
    Algorithm:
        for each BALL in BALLS:
            apply fn to BALL with args
    '''
    for ball in balls:
        fn(ball, *args)

# Examples:
fn_for_balls(balls0, set_pos, pos_x1, pos_y1)   # base case: empty BALLS, nothing happens
fn_for_balls(balls2, set_pos, pos_x0, pos_y0)   # every BALL moved to (0, 0)
fn_for_balls(balls2, set_vel, vel_x1, vel_y1)   # every BALL gets velocity (0.5, 0.5)
```

Use these sparingly. Most of the time, a single function applied to a single BALL should be called directly: `set_pos(ball, x, y)`. But using this wrapper is great for inserting evaluators and enforcing our invariants! Something that can (and does) bite after.

Sometimes we do need to access specific components of the ball outside of its scope. In these cases:

```
... ball.pos_x
... ball.pos_y
... ball.vel_x
... ball.vel_y
... ball.draw()
```

Notice this protocol also applies to dictionaries, only without the functions. In fact, the use of functions inside of dictionaries is a VIOLATION of our code practices.

### Template for dictionaries


```
# DD. USER
# user = {"NAME":str, "PHONE":str, "ON_DUTY":bool}
# interp. a user who works at a specific company
user0 = {"NAME":"Diego", "PHONE":"555-342-7445, "ON_DUTY":True}
user1 = {"NAME":"Raul", "PHONE":"555-362-1834, "ON_DUTY":False}


# Template for USER
... (user)
... user["NAME"]
... user["PHONE"]
... user["ON_DUTY"]

```

### Template for list and array-type objects

For array-type objects, we can use an iterator, whether it's a generator, for loop, while loop, etc.

```
Template for BALLS
for ball in balls:
    ... ball

```

### compounding templates

Building a template precedes the code. Let's take a look at how we can update the position of each BALL of BALLS found inside a GAME of GAMES using templates (my favorite part):

```
for game in games:
    for ball in game.balls:
        ball.simple_bounce()

# or

def update_balls(game):
    for ball in game.balls:
        ball.simple_bounce()

for game in games:
    update_balls(game)
```

However, when the proper functions are put in place, the use of templates shoudl cascade downstream as we place each template in their own helper function:

```
for game in games:
    game.update_balls()

# which can also be written as 
[game.update_balls() for game in games]
```




### docstrings and well commented code

When possible, function docstrings should have a small description of the algorithm to help the user understand what the implementation does. Such algorithmic disclosure can be presented as pseudocode. In addition, if they consume clearly defined data definitions, a small note in above their definition should declare the signature. An example:

```
# FD. shoot_ball
# Signature: COUNT, BALL -> COUNT
def shoot_ball(count:int) -> int:
    '''Instantiates a new ball at coordinates x,y for the ball closest to the player while adding a random velocity in x and y to the direction of the new ball
       Algorithm:
        Take BALL.x and BALL.y
        Instantiate a new BALL(BALL.x, BALL.y, random VEL_X, random VEL_Y)
        Increase COUNT by 1
        Return COUNT
    '''
    ...
```

The detailed structure of the data definitions encourages the implementation of clean and concise code.

## An introduction to the information domain

Each **atomic algorithm** represents one action, which takes input in and implements one clear action. Actions like:
- mapping one form of data into another
- filtering by selectively choosing elements by a criterion (including matching, slicing and deduplication)
- counting
- evaluating cases to find extremes (testing)
- Sorting (includes reversing, shuffling)
- appending (including concatenating)

A **compound algorithm** is the one that uses atomic algorithms to build its own behavior at a higher level of organization, delegating the protocols that are specific to atomic units to such elemental algorithms. Therefore, an algorithm that is aimed to sort and deduplicate should do so in two steps: first sort, then remove duplicates. 

If the separation of concerns results in performance loss, then this is a problem with the interpretation on the semantics of the problem, at which point we will **need to discuss together**. This will yield the highest rewards in the pursue of our goal.

Every code file (e.g. engine.py) will be assumed to be a compound algorithm, with clear data definitions and atomic algorithms inside. A .md file will always be accompanied by the .py file/compound algorithm, which will give a light introduction to the problem the algorithm is trying to solve, the semantics, a clear disclosure on the data definitions and invariants of each definition, and a clear, non technical explanation on what the algorithms do. After defining the problem to be solved and the declaration on the data definitions, we will define the atomic algorithms in plain english and finally do the same with the compound algorithm.

As follows:

```
module name: 
purpose: 

Introduction (up to 500 words). Explain what the problem in the information domain is so the user can read it and understand the role of this module at higher levels. Frame it in the context of the program's goal, then explain the details on what makes the module relevant.

Data definitions. Define the type of data, what it means, what it does. How does it connect to the problem being solved and how is it used by the compound algorithm?

Atomic algorithms: Title, what it does, what input does it take. What outputs/actions does it give? Easy to understand for a user reading through documents and code.

Compound algorithms: As in atomic algorithms, but clearly discloses how the atomic algorithms stack to form such compound. A diagram would be nice here, in addition to the paragraph.


```