(define (domain westbrick)
(:requirements :strips :equality)
(:predicates (clear ?obj)
             (objectat ?loc)
             (handempty)
             (holding ?obj)
             (robotat ?loc))

(:action pick
    :parameters (?obj ?loc)
    :precondition (and (objectat ?loc) (handempty) (clear ?obj))
    :effect (and  (holding ?obj) 
                  (not (clear ?obj)) (not (handempty)) (not (objectat ?loc))))

(:action place
    :parameters (?obj ?loc)
    :precondition (and (holding ?obj) (not (objectat ?loc)))
    :effect (and (handempty) (clear ?obj) (objectat ?loc)
                 (not (holding ?obj)))
)

(:action move
    :parameters (?loc)
    :precondition (and (not (robotat ?loc)))
    :effect (and (robotat ?loc))
)


)