from __future__ import annotations

from rest_framework import serializers

MAX_QUESTION_LENGTH = 2_000
MAX_TURN_LENGTH = 4_000
MAX_HISTORY_TURNS = 20


class TurnSerializer(serializers.Serializer):
    role = serializers.ChoiceField(choices=("user", "assistant"))
    content = serializers.CharField(max_length=MAX_TURN_LENGTH, allow_blank=False)


class TutorRequestSerializer(serializers.Serializer):
    question = serializers.CharField(max_length=MAX_QUESTION_LENGTH, allow_blank=False)
    history = TurnSerializer(many=True, required=False, default=list)

    def validate_history(self, value):
        if len(value) > MAX_HISTORY_TURNS:
            raise serializers.ValidationError(
                f"Only the last {MAX_HISTORY_TURNS} turns are used — trim the conversation."
            )
        return value
